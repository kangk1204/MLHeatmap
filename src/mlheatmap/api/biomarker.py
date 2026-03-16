"""Biomarker discovery API endpoints."""

from __future__ import annotations

import asyncio
import json
import traceback
from typing import Any

import numpy as np
from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse

from mlheatmap.api.session import SessionCancelledError
from mlheatmap.api.validation import acquire_session_lease_or_error
from mlheatmap.core.capabilities import MODEL_SPECS, get_model_capability, normalize_model_name

router = APIRouter(tags=["biomarker"])

DEG_TABLE_LIMIT = 5000
DEG_PLOT_TRUNCATION_THRESHOLD = 20000


def _stale_inputs_response():
    return JSONResponse(
        {"error": "Analysis inputs changed during computation. Rerun with the current settings."},
        status_code=409,
    )


def _cancelled_response(message: str = "Analysis cancelled at user request."):
    return JSONResponse({"error": message}, status_code=409)


def _build_sample_groups_from_parts(
    sample_names: list[str],
    groups: dict[str, list[str]],
    excluded_samples: list[str],
) -> tuple[dict[str, list[int]] | None, str | None]:
    if len(groups) < 2:
        return None, "Need at least 2 groups"

    excluded = set(excluded_samples)
    name_to_idx = {name: idx for idx, name in enumerate(sample_names)}
    sample_groups: dict[str, list[int]] = {}
    for group_name, group_samples in groups.items():
        indices = [
            name_to_idx[sample]
            for sample in group_samples
            if sample in name_to_idx and sample not in excluded
        ]
        if indices:
            sample_groups[group_name] = indices

    if len(sample_groups) < 2:
        return None, f"Only {len(sample_groups)} group(s) have valid samples after exclusion; need at least 2"

    for group_name, indices in sample_groups.items():
        if len(indices) < 2:
            return None, f"Group '{group_name}' has only {len(indices)} sample(s) after exclusion; need at least 2"

    return sample_groups, None


def _snapshot_biomarker_inputs(session) -> tuple[dict[str, Any] | None, JSONResponse | None]:
    with session.state_lock:
        if session.normalized is None:
            return None, JSONResponse({"error": "No normalized data"}, status_code=404)
        sample_groups, validation_error = _build_sample_groups_from_parts(
            list(session.sample_names),
            dict(session.groups),
            list(session.excluded_samples),
        )
        if validation_error:
            return None, JSONResponse({"error": validation_error}, status_code=400)

        return {
            "analysis_revision": session.analysis_revision,
            "expression": session.normalized,
            "gene_names": list(session.gene_names),
            "sample_groups": sample_groups,
        }, None


def _snapshot_deg_inputs(
    session,
    *,
    reference_group: str | None,
) -> tuple[dict[str, Any] | None, JSONResponse | None]:
    with session.state_lock:
        if session.normalized is None:
            return None, JSONResponse({"error": "No normalized data"}, status_code=404)
        if len(session.groups) != 2:
            return None, JSONResponse({"error": "DEG requires exactly 2 groups"}, status_code=400)
        excluded = set(session.excluded_samples)
        name_to_idx = {name: idx for idx, name in enumerate(session.sample_names)}
        sample_groups: dict[str, list[int]] = {}
        for group_name, group_samples in session.groups.items():
            indices = [
                name_to_idx[sample]
                for sample in group_samples
                if sample in name_to_idx and sample not in excluded
            ]
            if indices:
                sample_groups[group_name] = indices

        if len(sample_groups) != 2:
            missing = [group_name for group_name in session.groups if group_name not in sample_groups]
            if missing:
                return None, JSONResponse(
                    {"error": f"Group(s) {missing} have no valid samples after exclusion"},
                    status_code=400,
                )
            return None, JSONResponse(
                {"error": f"Only {len(sample_groups)} group(s) have valid samples after exclusion; need exactly 2"},
                status_code=400,
            )

        for group_name, indices in sample_groups.items():
            if len(indices) < 2:
                return None, JSONResponse(
                    {"error": f"Group '{group_name}' has only {len(indices)} sample(s) after exclusion; need at least 2"},
                    status_code=400,
                )

        if reference_group:
            if reference_group not in sample_groups:
                valid = ", ".join(sample_groups.keys())
                return None, JSONResponse(
                    {"error": f"Unknown reference group '{reference_group}'. Choose one of: {valid}"},
                    status_code=400,
                )
            comparison_group = next(group for group in sample_groups if group != reference_group)
            sample_groups = {
                comparison_group: sample_groups[comparison_group],
                reference_group: sample_groups[reference_group],
            }

        return {
            "analysis_revision": session.analysis_revision,
            "expression": session.normalized,
            "gene_names": list(session.gene_names),
            "sample_groups": sample_groups,
            "effect_size_data": session.deg_abundance,
            "effect_size_basis": session.deg_effect_size_basis,
            "normalization_method": session.norm_method,
        }, None


def _build_deg_response(result: dict[str, Any], normalization_method: str) -> dict[str, Any]:
    plot_results = result["results"]
    table_results = plot_results
    results_truncated = False

    if len(plot_results) > DEG_PLOT_TRUNCATION_THRESHOLD:
        table_results = plot_results[:DEG_TABLE_LIMIT]
        significant = [row for row in plot_results if row["direction"] != "ns"]
        seen = {row["gene"] for row in table_results}
        for row in significant:
            if row["gene"] not in seen:
                table_results.append(row)
                seen.add(row["gene"])
        results_truncated = len(table_results) != len(plot_results)

    return {
        "results": plot_results,
        "plot_results": plot_results,
        "table_results": table_results,
        "results_truncated": results_truncated,
        "result_counts": {
            "total": len(plot_results),
            "plot": len(plot_results),
            "table": len(table_results),
        },
        "summary": result["summary"],
        "group_names": result["group_names"],
        "comparison_group": result["comparison_group"],
        "reference_group": result["reference_group"],
        "thresholds": result["thresholds"],
        "method": result["method"],
        "pvalue_type": result.get("pvalue_type", "fdr"),
        "effect_size_basis": result.get("effect_size_basis", ""),
        "normalization_method": normalization_method,
    }


@router.get("/biomarker/stream")
async def biomarker_stream(
    request: Request,
    session_id: str = Query(...),
    n_top_genes: int = Query(20, ge=5, le=200),
    n_estimators: int = Query(500, ge=50, le=2000),
    cv_folds: int = Query(5, ge=2, le=10),
    model: str = Query("rf"),
    panel_method: str = Query("forward"),
):
    """Run biomarker analysis with SSE progress streaming."""
    lease, error, normalized_session_id = acquire_session_lease_or_error(request, session_id)
    if error is not None:
        return error
    session = lease.session

    snapshot, snapshot_error = _snapshot_biomarker_inputs(session)
    if snapshot_error is not None:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return snapshot_error

    normalized_model = normalize_model_name(model)
    model_capability = get_model_capability(normalized_model)
    if not model_capability["known"]:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        valid_models = ", ".join(MODEL_SPECS.keys())
        return JSONResponse(
            {"error": f"Unknown model: {model}. Choose from: {valid_models}"},
            status_code=400,
        )
    if not model_capability["available"]:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return JSONResponse({"error": model_capability["unavailable_reason"]}, status_code=400)

    panel_method = (panel_method or "forward").strip().lower()
    valid_panels = {"forward", "lasso", "stability", "mrmr"}
    if panel_method not in valid_panels:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return JSONResponse(
            {"error": f"Unknown panel method: {panel_method}. Choose from: {', '.join(sorted(valid_panels))}"},
            status_code=400,
        )

    async def event_generator():
        progress_queue: asyncio.Queue[dict[str, object]] = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def progress_callback(step, pct, msg):
            loop.call_soon_threadsafe(
                progress_queue.put_nowait,
                {"step": step, "pct": pct, "msg": msg},
            )

        from mlheatmap.core.biomarker import run_biomarker_analysis

        future = loop.run_in_executor(
            None,
            lambda: run_biomarker_analysis(
                expression=snapshot["expression"],
                gene_names=snapshot["gene_names"],
                sample_groups=snapshot["sample_groups"],
                n_top_genes=n_top_genes,
                n_estimators=n_estimators,
                cv_folds=cv_folds,
                model=normalized_model,
                panel_method=panel_method,
                progress_callback=progress_callback,
                cancel_check=lease.cancel_event.is_set,
            ),
        )

        try:
            while not future.done():
                if await request.is_disconnected():
                    lease.cancel_event.set()
                    return
                try:
                    progress = await asyncio.wait_for(progress_queue.get(), timeout=0.3)
                    yield f"event: progress\ndata: {json.dumps(progress)}\n\n"
                except asyncio.TimeoutError:
                    continue

            while not progress_queue.empty():
                try:
                    progress = progress_queue.get_nowait()
                    yield f"event: progress\ndata: {json.dumps(progress)}\n\n"
                except asyncio.QueueEmpty:
                    break

            result = future.result()
            with session.state_lock:
                if session.analysis_revision != snapshot["analysis_revision"]:
                    yield (
                        "event: app_error\ndata: "
                        f"{json.dumps({'detail': 'Analysis inputs changed during execution. Rerun the analysis.'})}\n\n"
                    )
                    return

                session.biomarker_results = result
                session.metadata["biomarker"] = {
                    "model": normalized_model,
                    "model_label": result.get("model", normalized_model),
                    "panel_method": panel_method,
                    "n_top_genes": n_top_genes,
                    "n_estimators": n_estimators,
                    "cv_folds_requested": result.get("cv_folds_requested", cv_folds),
                    "cv_folds_used": result.get(
                        "cv_folds_used",
                        min(cv_folds, min(len(v) for v in snapshot["sample_groups"].values())),
                    ),
                    "roc_evaluation": result.get("roc_evaluation"),
                    "panel_evaluation": result.get("optimal_combo", {}).get("evaluation"),
                    "shap_fallback_used": result.get("shap_fallback_used", False),
                    "shap_fallback_folds": result.get("shap_fallback_folds", []),
                }

            yield f"event: complete\ndata: {json.dumps(result, default=_json_safe)}\n\n"
        except SessionCancelledError as exc:
            yield f"event: app_error\ndata: {json.dumps({'detail': str(exc)})}\n\n"
        except Exception as exc:
            traceback.print_exc()
            yield f"event: app_error\ndata: {json.dumps({'detail': str(exc)})}\n\n"
        finally:
            request.app.state.sessions.end_use(lease.session_id, lease.operation_id)

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@router.get("/biomarker/deg")
async def deg_analysis(
    request: Request,
    session_id: str = Query(...),
    method: str = Query("wilcoxon"),
    log2fc_threshold: float = Query(1.0, ge=0),
    pvalue_threshold: float = Query(0.05, ge=0, le=1),
    use_raw_pvalue: bool = Query(False),
    reference_group: str | None = Query(None),
):
    """Run DEG analysis between groups."""
    lease, error, normalized_session_id = acquire_session_lease_or_error(request, session_id)
    if error is not None:
        return error
    session = lease.session

    snapshot, snapshot_error = _snapshot_deg_inputs(session, reference_group=reference_group)
    if snapshot_error is not None:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return snapshot_error

    valid_methods = ("wilcoxon", "ttest")
    if method not in valid_methods:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return JSONResponse(
            {"error": f"Unknown DEG method '{method}'. Use one of: {', '.join(valid_methods)}"},
            status_code=400,
        )

    from mlheatmap.core.deg import compute_deg

    try:
        result = await asyncio.to_thread(
            compute_deg,
            expression=snapshot["expression"],
            gene_names=snapshot["gene_names"],
            sample_groups=snapshot["sample_groups"],
            method=method,
            log2fc_threshold=log2fc_threshold,
            pvalue_threshold=pvalue_threshold,
            use_raw_pvalue=use_raw_pvalue,
            effect_size_data=snapshot["effect_size_data"],
            effect_size_basis=snapshot["effect_size_basis"],
            cancel_check=lease.cancel_event.is_set,
        )

        with session.state_lock:
            if session.analysis_revision != snapshot["analysis_revision"]:
                return _stale_inputs_response()

            session.deg_results = result
            session.metadata["deg"] = {
                "method": method,
                "pvalue_type": result.get("pvalue_type", "fdr"),
                "log2fc_threshold": log2fc_threshold,
                "pvalue_threshold": pvalue_threshold,
                "reference_group": result.get("reference_group", ""),
                "comparison_group": result.get("comparison_group", ""),
                "effect_size_basis": result.get("effect_size_basis", snapshot["effect_size_basis"]),
                "normalization_method": snapshot["normalization_method"],
            }

        response = _build_deg_response(result, snapshot["normalization_method"])
        return JSONResponse(response, headers={"Content-Type": "application/json"})
    except SessionCancelledError as exc:
        return _cancelled_response(str(exc))
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)


def _json_safe(obj):
    """JSON serializer for numpy types."""
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
