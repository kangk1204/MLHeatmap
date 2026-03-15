"""Biomarker discovery API endpoints."""

from __future__ import annotations

import asyncio
import json
import traceback

import numpy as np
from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse

from mlheatmap.core.capabilities import MODEL_SPECS, get_model_capability, normalize_model_name

router = APIRouter(tags=["biomarker"])


def _stale_inputs_response():
    return JSONResponse(
        {"error": "Analysis inputs changed during computation. Rerun with the current settings."},
        status_code=409,
    )


def _build_sample_groups(session) -> tuple[dict[str, list[int]] | None, str | None]:
    if len(session.groups) < 2:
        return None, "Need at least 2 groups"

    sample_groups: dict[str, list[int]] = {}
    for group, samples in session.groups.items():
        for sample in samples:
            if sample in session.sample_names and sample not in session.excluded_samples:
                sample_groups.setdefault(group, []).append(session.sample_names.index(sample))

    if len(sample_groups) < 2:
        return None, f"Only {len(sample_groups)} group(s) have valid samples after exclusion; need at least 2"

    for group_name, indices in sample_groups.items():
        if len(indices) < 2:
            return None, f"Group '{group_name}' has only {len(indices)} sample(s) after exclusion; need at least 2"

    return sample_groups, None


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
    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    sample_groups, validation_error = _build_sample_groups(session)
    if validation_error:
        return JSONResponse({"error": validation_error}, status_code=400)

    normalized_model = normalize_model_name(model)
    model_capability = get_model_capability(normalized_model)
    if not model_capability["known"]:
        valid_models = ", ".join(MODEL_SPECS.keys())
        return JSONResponse(
            {"error": f"Unknown model: {model}. Choose from: {valid_models}"},
            status_code=400,
        )
    if not model_capability["available"]:
        return JSONResponse({"error": model_capability["unavailable_reason"]}, status_code=400)

    valid_panels = {"forward", "lasso", "stability", "mrmr"}
    if panel_method not in valid_panels:
        return JSONResponse(
            {"error": f"Unknown panel method: {panel_method}. Choose from: {', '.join(sorted(valid_panels))}"},
            status_code=400,
        )

    async def event_generator():
        start_revision = session.analysis_revision
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
                expression=session.normalized,
                gene_names=session.gene_names,
                sample_groups=sample_groups,
                n_top_genes=n_top_genes,
                n_estimators=n_estimators,
                cv_folds=cv_folds,
                model=normalized_model,
                panel_method=panel_method,
                progress_callback=progress_callback,
            ),
        )

        while not future.done():
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

        try:
            result = future.result()
            if session.analysis_revision != start_revision:
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
                "cv_folds_requested": cv_folds,
                "cv_folds_used": min(cv_folds, min(len(v) for v in sample_groups.values())),
                "roc_evaluation": result.get("roc_evaluation"),
                "panel_evaluation": result.get("optimal_combo", {}).get("evaluation"),
            }
            yield f"event: complete\ndata: {json.dumps(result, default=_json_safe)}\n\n"
        except Exception as exc:
            traceback.print_exc()
            yield f"event: app_error\ndata: {json.dumps({'detail': str(exc)})}\n\n"

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
    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)
    start_revision = session.analysis_revision

    if len(session.groups) != 2:
        return JSONResponse({"error": "DEG requires exactly 2 groups"}, status_code=400)

    sample_groups: dict[str, list[int]] = {}
    for group, samples in session.groups.items():
        for sample in samples:
            if sample in session.sample_names and sample not in session.excluded_samples:
                sample_groups.setdefault(group, []).append(session.sample_names.index(sample))

    if len(sample_groups) != 2:
        missing = [group for group in session.groups if group not in sample_groups]
        return JSONResponse(
            {"error": f"Group(s) {missing} have no valid samples after exclusion"},
            status_code=400,
        )

    for group_name, indices in sample_groups.items():
        if len(indices) < 2:
            return JSONResponse(
                {"error": f"Group '{group_name}' has only {len(indices)} sample(s) after exclusion; need at least 2"},
                status_code=400,
            )

    valid_methods = ("wilcoxon", "ttest")
    if method not in valid_methods:
        return JSONResponse(
            {"error": f"Unknown DEG method '{method}'. Use one of: {', '.join(valid_methods)}"},
            status_code=400,
        )

    if reference_group:
        if reference_group not in sample_groups:
            valid = ", ".join(sample_groups.keys())
            return JSONResponse(
                {"error": f"Unknown reference group '{reference_group}'. Choose one of: {valid}"},
                status_code=400,
            )
        comparison = next(group for group in sample_groups if group != reference_group)
        sample_groups = {
            comparison: sample_groups[comparison],
            reference_group: sample_groups[reference_group],
        }

    from mlheatmap.core.deg import compute_deg

    result = compute_deg(
        expression=session.normalized,
        gene_names=session.gene_names,
        sample_groups=sample_groups,
        method=method,
        log2fc_threshold=log2fc_threshold,
        pvalue_threshold=pvalue_threshold,
        use_raw_pvalue=use_raw_pvalue,
        effect_size_data=session.deg_abundance,
        effect_size_basis=session.deg_effect_size_basis,
    )

    if session.analysis_revision != start_revision:
        return _stale_inputs_response()

    session.deg_results = result
    session.metadata["deg"] = {
        "method": method,
        "pvalue_type": result.get("pvalue_type", "fdr"),
        "log2fc_threshold": log2fc_threshold,
        "pvalue_threshold": pvalue_threshold,
        "reference_group": result.get("reference_group", ""),
        "comparison_group": result.get("comparison_group", ""),
        "effect_size_basis": result.get("effect_size_basis", session.deg_effect_size_basis),
        "normalization_method": session.norm_method,
    }

    all_results = result["results"]
    if len(all_results) > 20000:
        top_results = all_results[:5000]
        significant = [row for row in all_results if row["direction"] != "ns"]
        seen = {row["gene"] for row in top_results}
        for row in significant:
            if row["gene"] not in seen:
                top_results.append(row)
        all_results = top_results

    return JSONResponse(
        {
            "results": all_results,
            "summary": result["summary"],
            "group_names": result["group_names"],
            "comparison_group": result["comparison_group"],
            "reference_group": result["reference_group"],
            "thresholds": result["thresholds"],
            "method": result["method"],
            "pvalue_type": result.get("pvalue_type", "fdr"),
            "effect_size_basis": result.get("effect_size_basis", session.deg_effect_size_basis),
            "normalization_method": session.norm_method,
        },
        headers={"Content-Type": "application/json"},
    )


def _json_safe(obj):
    """JSON serializer for numpy types."""
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
