"""Hierarchical clustering for heatmap."""

import logging
import numpy as np
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist

logger = logging.getLogger(__name__)


def select_top_variable_genes(
    expression: np.ndarray,
    gene_names: list[str],
    top_n: int,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Select top-N most variable genes by variance.

    Returns:
        Tuple of (expression_subset, gene_names_subset, indices)
    """
    n_genes = expression.shape[0]
    top_n = min(top_n, n_genes)
    variances = np.var(expression, axis=1)
    top_idx = np.argsort(variances)[-top_n:]
    expr_sub = expression[top_idx]
    genes_sub = [gene_names[i] for i in top_idx]
    return expr_sub, genes_sub, top_idx


def compute_heatmap_data(
    expression: np.ndarray,
    gene_names: list[str],
    sample_names: list[str],
    top_n: int = 500,
    distance: str = "correlation",
    method: str = "average",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
) -> dict:
    """Compute clustered heatmap data with dendrograms.

    Args:
        expression: Normalized matrix (genes x samples)
        gene_names: Gene names
        sample_names: Sample names
        top_n: Number of most variable genes to show
        distance: Distance metric
        method: Linkage method
        cluster_rows: Whether to cluster rows (genes)
        cluster_cols: Whether to cluster columns (samples)

    Returns:
        Dict with z-scored matrix, ordered names, dendrogram data
    """
    n_genes, n_samples = expression.shape

    # Replace non-finite values (from log2(0) = -inf, etc.)
    expression = np.nan_to_num(expression, nan=0.0, posinf=0.0, neginf=0.0)

    # Select top-N variable genes
    expr_sub, genes_sub, _ = select_top_variable_genes(expression, gene_names, top_n)
    top_n = len(genes_sub)

    # Row clustering (genes)
    row_dendro_data = None
    row_order = np.arange(top_n)
    if cluster_rows and top_n > 1:
        try:
            if distance == "correlation":
                row_dist = pdist(expr_sub, metric="correlation")
                row_dist = np.nan_to_num(row_dist, nan=0.0)
            else:
                row_dist = pdist(expr_sub, metric=distance)
            row_link = linkage(row_dist, method=method)
            row_d = dendrogram(row_link, no_plot=True)
            row_order = leaves_list(row_link)
            row_dendro_data = {
                "icoord": row_d["icoord"],
                "dcoord": row_d["dcoord"],
            }
        except Exception as e:
            logger.warning("Row clustering failed: %s", e)

    # Column clustering (samples)
    col_dendro_data = None
    col_order = np.arange(n_samples)
    if cluster_cols and n_samples > 1:
        try:
            col_metric = "correlation" if distance == "correlation" else distance
            col_dist = pdist(expr_sub.T, metric=col_metric)
            col_dist = np.nan_to_num(col_dist, nan=0.0)
            col_link = linkage(col_dist, method=method)
            col_d = dendrogram(col_link, no_plot=True)
            col_order = leaves_list(col_link)
            col_dendro_data = {
                "icoord": col_d["icoord"],
                "dcoord": col_d["dcoord"],
            }
        except Exception as e:
            logger.warning("Column clustering failed: %s", e)

    # Reorder
    ordered = expr_sub[row_order][:, col_order]

    # Z-score rows
    row_means = np.mean(ordered, axis=1, keepdims=True)
    row_stds = np.std(ordered, axis=1, keepdims=True)
    row_stds = np.maximum(row_stds, 1e-6)
    z_scored = (ordered - row_means) / row_stds

    # Replace any remaining NaN/inf and clip extreme values
    z_scored = np.nan_to_num(z_scored, nan=0.0, posinf=3.0, neginf=-3.0)
    z_scored = np.clip(z_scored, -3, 3)

    return {
        "z": z_scored.tolist(),
        "x": [sample_names[i] for i in col_order],
        "y": [genes_sub[i] for i in row_order],
        "row_dendrogram": row_dendro_data,
        "col_dendrogram": col_dendro_data,
        "n_total_genes": n_genes,
        "n_shown_genes": top_n,
    }
