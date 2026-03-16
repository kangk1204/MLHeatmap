"""Heatmap API endpoints."""

from __future__ import annotations

import asyncio
import copy

import numpy as np
from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, Response

from mlheatmap.api.session import SessionCancelledError
from mlheatmap.api.validation import acquire_session_lease_or_error
from mlheatmap.core.clustering import validate_heatmap_params

router = APIRouter(tags=["heatmap"])

VALID_IMAGE_FORMATS = {"png", "svg"}


def _stale_inputs_response():
    return JSONResponse(
        {"error": "Analysis inputs changed during computation. Rerun with the current settings."},
        status_code=409,
    )


def _cancelled_response(message: str = "Analysis cancelled at user request."):
    return JSONResponse({"error": message}, status_code=409)


def _validate_request_params(distance: str, linkage: str, *, fmt: str | None = None) -> JSONResponse | None:
    try:
        validate_heatmap_params(distance, linkage)
    except ValueError as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)

    if fmt is not None and fmt not in VALID_IMAGE_FORMATS:
        valid = ", ".join(sorted(VALID_IMAGE_FORMATS))
        return JSONResponse({"error": f"Unknown format '{fmt}'. Choose one of: {valid}"}, status_code=400)

    return None


def _order_samples_by_groups(
    expression: np.ndarray,
    sample_names: list[str],
    groups: dict[str, list[str]],
    cluster_cols: bool,
) -> tuple[np.ndarray, list[str]]:
    """Keep samples grouped when column clustering is disabled."""
    if cluster_cols or expression.shape[1] <= 1 or not groups:
        return expression, sample_names

    name_to_idx = {name: idx for idx, name in enumerate(sample_names)}
    ordered_indices = []
    seen = set()

    for group_samples in groups.values():
        for sample in group_samples:
            idx = name_to_idx.get(sample)
            if idx is not None and idx not in seen:
                ordered_indices.append(idx)
                seen.add(idx)

    ordered_indices.extend(idx for idx in range(len(sample_names)) if idx not in seen)
    return expression[:, ordered_indices], [sample_names[idx] for idx in ordered_indices]


def _snapshot_heatmap_inputs(
    session,
    *,
    require_biomarker: bool = False,
    require_deg: bool = False,
) -> tuple[dict | None, JSONResponse | None]:
    with session.state_lock:
        if session.normalized is None:
            return None, JSONResponse({"error": "No normalized data"}, status_code=404)
        if require_biomarker and not session.biomarker_results:
            return None, JSONResponse({"error": "Run biomarker analysis first"}, status_code=400)
        if require_deg and not session.deg_results:
            return None, JSONResponse({"error": "Run DEG analysis first"}, status_code=400)

        return {
            "analysis_revision": session.analysis_revision,
            "expression": np.array(session.normalized, copy=True),
            "gene_names": list(session.gene_names),
            "sample_names": list(session.sample_names),
            "groups": copy.deepcopy(session.groups),
            "excluded_samples": list(session.excluded_samples),
            "biomarker_results": copy.deepcopy(session.biomarker_results) if require_biomarker else None,
            "deg_results": copy.deepcopy(session.deg_results) if require_deg else None,
        }, None


def _filtered_expression(
    snapshot: dict,
    cluster_cols: bool,
) -> tuple[np.ndarray | None, list[str] | None, JSONResponse | None]:
    sample_mask = [
        idx for idx, sample in enumerate(snapshot["sample_names"]) if sample not in snapshot["excluded_samples"]
    ]
    if not sample_mask:
        return None, None, JSONResponse(
            {"error": "No samples remain after exclusion. Re-include at least one sample."},
            status_code=400,
        )
    expression = snapshot["expression"][:, sample_mask]
    sample_names = [snapshot["sample_names"][idx] for idx in sample_mask]
    expression, sample_names = _order_samples_by_groups(expression, sample_names, snapshot["groups"], cluster_cols)
    if expression.shape[1] == 0:
        return None, None, JSONResponse(
            {"error": "No samples remain after exclusion. Re-include at least one sample."},
            status_code=400,
        )
    return expression, sample_names, None


@router.get("/heatmap")
async def get_heatmap(
    request: Request,
    session_id: str = Query(...),
    top_n: int = Query(500, ge=10, le=10000),
    distance: str = Query("correlation"),
    linkage: str = Query("average"),
    color_scale: str = Query("RdBu_r"),
    cluster_rows: bool = Query(True),
    cluster_cols: bool = Query(True),
):
    """Compute clustered heatmap data."""
    from mlheatmap.core.clustering import compute_heatmap_data

    error = _validate_request_params(distance, linkage)
    if error is not None:
        return error

    lease, acquire_error, normalized_session_id = acquire_session_lease_or_error(request, session_id)
    if acquire_error is not None:
        return acquire_error
    session = lease.session

    snapshot, snapshot_error = _snapshot_heatmap_inputs(session)
    if snapshot_error is not None:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return snapshot_error

    try:
        expression, sample_names, filtered_error = _filtered_expression(snapshot, cluster_cols)
        if filtered_error is not None:
            return filtered_error
        result = await asyncio.to_thread(
            compute_heatmap_data,
            expression=expression,
            gene_names=snapshot["gene_names"],
            sample_names=sample_names,
            top_n=min(top_n, expression.shape[0]),
            distance=distance,
            method=linkage,
            cluster_rows=cluster_rows,
            cluster_cols=cluster_cols,
            cancel_check=lease.cancel_event.is_set,
        )

        with session.state_lock:
            if session.analysis_revision != snapshot["analysis_revision"]:
                return _stale_inputs_response()
            session.heatmap_data = result
            session.heatmap_color_scale = color_scale
            groups = copy.deepcopy(session.groups)

        response = dict(result)
        response["groups"] = groups
        response["color_scale"] = color_scale
        return response
    except SessionCancelledError as exc:
        return _cancelled_response(str(exc))
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)


@router.get("/heatmap/render")
async def render_heatmap(
    request: Request,
    session_id: str = Query(...),
    top_n: int = Query(500, ge=10, le=60000),
    distance: str = Query("correlation"),
    linkage: str = Query("average"),
    color_scale: str = Query("RdBu_r"),
    cluster_rows: bool = Query(True),
    cluster_cols: bool = Query(True),
    fmt: str = Query("png"),
    dpi: int = Query(150, ge=72, le=600),
):
    """Server-side rendered heatmap image for large gene sets."""
    from mlheatmap.core.server_render import render_heatmap_image

    error = _validate_request_params(distance, linkage, fmt=fmt)
    if error is not None:
        return error

    lease, acquire_error, normalized_session_id = acquire_session_lease_or_error(request, session_id)
    if acquire_error is not None:
        return acquire_error
    session = lease.session

    snapshot, snapshot_error = _snapshot_heatmap_inputs(session)
    if snapshot_error is not None:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return snapshot_error

    try:
        expression, sample_names, filtered_error = _filtered_expression(snapshot, cluster_cols)
        if filtered_error is not None:
            return filtered_error
        image_bytes = await asyncio.to_thread(
            render_heatmap_image,
            expression=expression,
            gene_names=snapshot["gene_names"],
            sample_names=sample_names,
            groups=snapshot["groups"],
            top_n=min(top_n, expression.shape[0]),
            distance=distance,
            method=linkage,
            color_scale=color_scale,
            cluster_rows=cluster_rows,
            cluster_cols=cluster_cols,
            fmt=fmt,
            dpi=dpi,
            cancel_check=lease.cancel_event.is_set,
        )

        with session.state_lock:
            if session.analysis_revision != snapshot["analysis_revision"]:
                return _stale_inputs_response()
            session.heatmap_color_scale = color_scale

        media_type = "image/svg+xml" if fmt == "svg" else "image/png"
        return Response(image_bytes, media_type=media_type, headers={"Cache-Control": "no-store"})
    except SessionCancelledError as exc:
        return _cancelled_response(str(exc))
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)


@router.get("/heatmap/shap")
async def get_shap_heatmap(
    request: Request,
    session_id: str = Query(...),
    top_n: int = Query(20, ge=5, le=200),
    distance: str = Query("correlation"),
    linkage: str = Query("average"),
    color_scale: str = Query("RdBu_r"),
    cluster_rows: bool = Query(True),
    cluster_cols: bool = Query(True),
):
    """Compute heatmap using SHAP-ranked top genes from biomarker results."""
    from mlheatmap.core.clustering import compute_heatmap_data

    error = _validate_request_params(distance, linkage)
    if error is not None:
        return error

    lease, acquire_error, normalized_session_id = acquire_session_lease_or_error(request, session_id)
    if acquire_error is not None:
        return acquire_error
    session = lease.session

    snapshot, snapshot_error = _snapshot_heatmap_inputs(session, require_biomarker=True)
    if snapshot_error is not None:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return snapshot_error

    biomarker_results = snapshot["biomarker_results"]
    top_genes = biomarker_results["top_genes"][:top_n]
    shap_gene_symbols = [gene["symbol"] for gene in top_genes]
    shap_values = [gene["shap_mean_abs"] for gene in top_genes]

    gene_name_to_idx = {name: idx for idx, name in enumerate(snapshot["gene_names"])}
    gene_indices = []
    found_symbols = []
    found_shap = []
    for symbol, shap_value in zip(shap_gene_symbols, shap_values):
        idx = gene_name_to_idx.get(symbol)
        if idx is not None:
            gene_indices.append(idx)
            found_symbols.append(symbol)
            found_shap.append(shap_value)

    if not gene_indices:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return JSONResponse({"error": "No matching genes found"}, status_code=400)

    try:
        expression, sample_names, filtered_error = _filtered_expression(snapshot, cluster_cols)
        if filtered_error is not None:
            return filtered_error
        subset_expression = expression[np.array(gene_indices)]
        result = await asyncio.to_thread(
            compute_heatmap_data,
            expression=subset_expression,
            gene_names=found_symbols,
            sample_names=sample_names,
            top_n=len(found_symbols),
            distance=distance,
            method=linkage,
            cluster_rows=cluster_rows,
            cluster_cols=cluster_cols,
            cancel_check=lease.cancel_event.is_set,
        )

        with session.state_lock:
            if session.analysis_revision != snapshot["analysis_revision"]:
                return _stale_inputs_response()
            session.heatmap_color_scale = color_scale
            groups = copy.deepcopy(session.groups)

        response = dict(result)
        response["groups"] = groups
        response["color_scale"] = color_scale
        response["model"] = biomarker_results.get("model", "Random Forest")
        shap_map = dict(zip(found_symbols, found_shap))
        response["shap_values"] = [shap_map.get(gene, 0) for gene in result["y"]]
        return response
    except SessionCancelledError as exc:
        return _cancelled_response(str(exc))
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)


@router.get("/heatmap/deg")
async def get_deg_heatmap(
    request: Request,
    session_id: str = Query(...),
    top_n: int = Query(30, ge=5, le=5000),
    distance: str = Query("correlation"),
    linkage: str = Query("average"),
    color_scale: str = Query("RdBu_r"),
    cluster_rows: bool = Query(True),
    cluster_cols: bool = Query(True),
):
    """Compute heatmap using DEG-ranked top genes."""
    from mlheatmap.core.clustering import compute_heatmap_data

    error = _validate_request_params(distance, linkage)
    if error is not None:
        return error

    lease, acquire_error, normalized_session_id = acquire_session_lease_or_error(request, session_id)
    if acquire_error is not None:
        return acquire_error
    session = lease.session

    snapshot, snapshot_error = _snapshot_heatmap_inputs(session, require_deg=True)
    if snapshot_error is not None:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return snapshot_error

    deg_results = snapshot["deg_results"]
    deg_genes = deg_results["results"][:top_n]
    deg_symbols = [gene["gene"] for gene in deg_genes]
    deg_log2fc = [gene["log2fc"] for gene in deg_genes]
    deg_neglog10p = [gene["neg_log10_p"] for gene in deg_genes]

    gene_name_to_idx = {name: idx for idx, name in enumerate(snapshot["gene_names"])}
    gene_indices = []
    found_symbols = []
    found_log2fc = []
    found_neglog10p = []
    for symbol, log2fc, neglog10p in zip(deg_symbols, deg_log2fc, deg_neglog10p):
        idx = gene_name_to_idx.get(symbol)
        if idx is not None:
            gene_indices.append(idx)
            found_symbols.append(symbol)
            found_log2fc.append(log2fc)
            found_neglog10p.append(neglog10p)

    if not gene_indices:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
        return JSONResponse({"error": "No matching genes found"}, status_code=400)

    try:
        expression, sample_names, filtered_error = _filtered_expression(snapshot, cluster_cols)
        if filtered_error is not None:
            return filtered_error
        subset_expression = expression[np.array(gene_indices)]
        result = await asyncio.to_thread(
            compute_heatmap_data,
            expression=subset_expression,
            gene_names=found_symbols,
            sample_names=sample_names,
            top_n=len(found_symbols),
            distance=distance,
            method=linkage,
            cluster_rows=cluster_rows,
            cluster_cols=cluster_cols,
            cancel_check=lease.cancel_event.is_set,
        )

        with session.state_lock:
            if session.analysis_revision != snapshot["analysis_revision"]:
                return _stale_inputs_response()
            session.heatmap_color_scale = color_scale
            groups = copy.deepcopy(session.groups)

        response = dict(result)
        response["groups"] = groups
        response["color_scale"] = color_scale
        response["method"] = deg_results.get("method", "wilcoxon")
        response["comparison_group"] = deg_results.get("comparison_group", "")
        response["reference_group"] = deg_results.get("reference_group", "")
        response["effect_size_basis"] = deg_results.get("effect_size_basis", "")

        fc_map = dict(zip(found_symbols, found_log2fc))
        nlp_map = dict(zip(found_symbols, found_neglog10p))
        response["log2fc_values"] = [fc_map.get(gene, 0) for gene in result["y"]]
        response["neglog10p_values"] = [nlp_map.get(gene, 0) for gene in result["y"]]
        return response
    except SessionCancelledError as exc:
        return _cancelled_response(str(exc))
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
