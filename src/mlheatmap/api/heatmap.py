"""Heatmap API endpoints."""

import asyncio
import numpy as np

from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, Response

router = APIRouter(tags=["heatmap"])


def _order_samples_by_groups(
    expression: np.ndarray,
    sample_names: list[str],
    groups: dict,
    cluster_cols: bool,
) -> tuple[np.ndarray, list[str]]:
    """Keep samples grouped when column clustering is disabled.

    The group order follows the current group-definition order from the UI.
    Samples not assigned to any group are appended at the end in their
    original order.
    """
    if cluster_cols or expression.shape[1] <= 1 or not groups:
        return expression, sample_names

    name_to_idx = {name: i for i, name in enumerate(sample_names)}
    ordered_indices = []
    seen = set()

    for group_samples in groups.values():
        for sample in group_samples:
            idx = name_to_idx.get(sample)
            if idx is not None and idx not in seen:
                ordered_indices.append(idx)
                seen.add(idx)

    ordered_indices.extend(i for i in range(len(sample_names)) if i not in seen)

    return expression[:, ordered_indices], [sample_names[i] for i in ordered_indices]


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

    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    # Exclude samples if any
    sample_mask = [
        i for i, s in enumerate(session.sample_names)
        if s not in session.excluded_samples
    ]
    expression = session.normalized[:, sample_mask]
    sample_names = [session.sample_names[i] for i in sample_mask]
    expression, sample_names = _order_samples_by_groups(
        expression, sample_names, session.groups, cluster_cols
    )

    result = compute_heatmap_data(
        expression=expression,
        gene_names=session.gene_names,
        sample_names=sample_names,
        top_n=min(top_n, expression.shape[0]),
        distance=distance,
        method=linkage,
        cluster_rows=cluster_rows,
        cluster_cols=cluster_cols,
    )

    session.heatmap_data = result
    session.heatmap_color_scale = color_scale

    # Return a copy with groups/color_scale to avoid mutating stored data
    response = dict(result)
    response["groups"] = session.groups
    response["color_scale"] = color_scale

    return response


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

    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    # Exclude samples
    sample_mask = [
        i for i, s in enumerate(session.sample_names)
        if s not in session.excluded_samples
    ]
    expression = session.normalized[:, sample_mask]
    sample_names = [session.sample_names[i] for i in sample_mask]
    expression, sample_names = _order_samples_by_groups(
        expression, sample_names, session.groups, cluster_cols
    )

    session.heatmap_color_scale = color_scale

    try:
        image_bytes = await asyncio.to_thread(
            render_heatmap_image,
            expression=expression,
            gene_names=session.gene_names,
            sample_names=sample_names,
            groups=session.groups,
            top_n=min(top_n, expression.shape[0]),
            distance=distance,
            method=linkage,
            color_scale=color_scale,
            cluster_rows=cluster_rows,
            cluster_cols=cluster_cols,
            fmt=fmt,
            dpi=dpi,
        )
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=500)

    media = "image/svg+xml" if fmt == "svg" else "image/png"
    return Response(image_bytes, media_type=media, headers={"Cache-Control": "no-store"})


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

    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    if not session.biomarker_results:
        return JSONResponse({"error": "Run biomarker analysis first"}, status_code=400)

    # Get SHAP-ranked gene symbols
    top_genes = session.biomarker_results["top_genes"][:top_n]
    shap_gene_symbols = [g["symbol"] for g in top_genes]
    shap_values = [g["shap_mean_abs"] for g in top_genes]

    # Find indices of these genes in the expression matrix
    gene_name_to_idx = {name: i for i, name in enumerate(session.gene_names)}
    gene_indices = []
    found_symbols = []
    found_shap = []
    for sym, shap_val in zip(shap_gene_symbols, shap_values):
        if sym in gene_name_to_idx:
            gene_indices.append(gene_name_to_idx[sym])
            found_symbols.append(sym)
            found_shap.append(shap_val)

    if len(gene_indices) == 0:
        return JSONResponse({"error": "No matching genes found"}, status_code=400)

    # Exclude samples
    sample_mask = [
        i for i, s in enumerate(session.sample_names)
        if s not in session.excluded_samples
    ]
    expression = session.normalized[np.array(gene_indices)][:, sample_mask]
    sample_names = [session.sample_names[i] for i in sample_mask]
    expression, sample_names = _order_samples_by_groups(
        expression, sample_names, session.groups, cluster_cols
    )

    result = compute_heatmap_data(
        expression=expression,
        gene_names=found_symbols,
        sample_names=sample_names,
        top_n=len(found_symbols),  # use all (already filtered)
        distance=distance,
        method=linkage,
        cluster_rows=cluster_rows,
        cluster_cols=cluster_cols,
    )

    session.heatmap_color_scale = color_scale

    response = dict(result)
    response["groups"] = session.groups
    response["color_scale"] = color_scale
    response["model"] = session.biomarker_results.get("model", "Random Forest")

    # Reorder SHAP values to match clustered y order
    shap_map = dict(zip(found_symbols, found_shap))
    response["shap_values"] = [shap_map.get(gene, 0) for gene in result["y"]]

    return response


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

    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    if not hasattr(session, "deg_results") or not session.deg_results:
        return JSONResponse({"error": "Run DEG analysis first"}, status_code=400)

    # Get DEG-ranked gene symbols (sorted by adj_pvalue)
    deg_genes = session.deg_results["results"][:top_n]
    deg_symbols = [g["gene"] for g in deg_genes]
    deg_log2fc = [g["log2fc"] for g in deg_genes]
    deg_neglog10p = [g["neg_log10_p"] for g in deg_genes]

    # Find indices in expression matrix
    gene_name_to_idx = {name: i for i, name in enumerate(session.gene_names)}
    gene_indices = []
    found_symbols = []
    found_log2fc = []
    found_neglog10p = []
    for sym, fc, nlp in zip(deg_symbols, deg_log2fc, deg_neglog10p):
        if sym in gene_name_to_idx:
            gene_indices.append(gene_name_to_idx[sym])
            found_symbols.append(sym)
            found_log2fc.append(fc)
            found_neglog10p.append(nlp)

    if len(gene_indices) == 0:
        return JSONResponse({"error": "No matching genes found"}, status_code=400)

    # Exclude samples
    sample_mask = [
        i for i, s in enumerate(session.sample_names)
        if s not in session.excluded_samples
    ]
    expression = session.normalized[np.array(gene_indices)][:, sample_mask]
    sample_names = [session.sample_names[i] for i in sample_mask]
    expression, sample_names = _order_samples_by_groups(
        expression, sample_names, session.groups, cluster_cols
    )

    result = compute_heatmap_data(
        expression=expression,
        gene_names=found_symbols,
        sample_names=sample_names,
        top_n=len(found_symbols),
        distance=distance,
        method=linkage,
        cluster_rows=cluster_rows,
        cluster_cols=cluster_cols,
    )

    session.heatmap_color_scale = color_scale

    response = dict(result)
    response["groups"] = session.groups
    response["color_scale"] = color_scale
    response["method"] = session.deg_results.get("method", "wilcoxon")
    response["comparison_group"] = session.deg_results.get("comparison_group", "")
    response["reference_group"] = session.deg_results.get("reference_group", "")

    # Reorder DEG values to match clustered y order
    fc_map = dict(zip(found_symbols, found_log2fc))
    nlp_map = dict(zip(found_symbols, found_neglog10p))
    response["log2fc_values"] = [fc_map.get(gene, 0) for gene in result["y"]]
    response["neglog10p_values"] = [nlp_map.get(gene, 0) for gene in result["y"]]

    return response
