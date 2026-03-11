"""Heatmap API endpoints."""

from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse

router = APIRouter(tags=["heatmap"])


@router.get("/heatmap")
async def get_heatmap(
    request: Request,
    session_id: str = Query(...),
    top_n: int = Query(500, ge=10, le=5000),
    distance: str = Query("correlation"),
    linkage: str = Query("average"),
    color_scale: str = Query("RdBu_r"),
    cluster_rows: bool = Query(True),
    cluster_cols: bool = Query(True),
):
    """Compute clustered heatmap data."""
    import numpy as np
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

    # Return a copy with groups/color_scale to avoid mutating stored data
    response = dict(result)
    response["groups"] = session.groups
    response["color_scale"] = color_scale

    return response


@router.get("/heatmap/shap")
async def get_shap_heatmap(
    request: Request,
    session_id: str = Query(...),
    top_n: int = Query(20, ge=5, le=100),
    distance: str = Query("correlation"),
    linkage: str = Query("average"),
    color_scale: str = Query("RdBu_r"),
    cluster_rows: bool = Query(True),
    cluster_cols: bool = Query(True),
):
    """Compute heatmap using SHAP-ranked top genes from biomarker results."""
    import numpy as np
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

    response = dict(result)
    response["groups"] = session.groups
    response["color_scale"] = color_scale
    response["model"] = session.biomarker_results.get("model", "Random Forest")

    # Reorder SHAP values to match clustered y order
    shap_map = dict(zip(found_symbols, found_shap))
    response["shap_values"] = [shap_map.get(gene, 0) for gene in result["y"]]

    return response
