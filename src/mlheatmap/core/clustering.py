"""Hierarchical clustering for heatmap."""

from __future__ import annotations

import logging

import numpy as np
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist

from mlheatmap.core.cancellation import raise_if_cancelled

logger = logging.getLogger(__name__)

VALID_DISTANCE_METRICS = {"correlation", "euclidean", "cityblock", "cosine"}
VALID_LINKAGE_METHODS = {"average", "complete", "single", "ward"}
ZSCORE_CLIP_LIMIT = 3.0


def validate_heatmap_params(distance: str, method: str) -> None:
    """Validate clustering parameters before expensive work begins."""
    if distance not in VALID_DISTANCE_METRICS:
        valid = ", ".join(sorted(VALID_DISTANCE_METRICS))
        raise ValueError(f"Unknown distance metric '{distance}'. Choose one of: {valid}")
    if method not in VALID_LINKAGE_METHODS:
        valid = ", ".join(sorted(VALID_LINKAGE_METHODS))
        raise ValueError(f"Unknown linkage method '{method}'. Choose one of: {valid}")
    if method == "ward" and distance != "euclidean":
        raise ValueError("Ward linkage requires euclidean distance.")


def select_top_variable_genes(
    expression: np.ndarray,
    gene_names: list[str],
    top_n: int,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Select top-N most variable genes by variance."""
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
    cancel_check=None,
) -> dict:
    """Compute clustered heatmap data with dendrograms."""
    validate_heatmap_params(distance, method)
    if expression.ndim != 2:
        raise ValueError("Expression matrix must be 2-dimensional.")
    if expression.shape[0] == 0:
        raise ValueError("No genes are available for heatmap rendering.")
    if expression.shape[1] == 0 or not sample_names:
        raise ValueError("No samples remain after exclusion. Re-include at least one sample.")

    n_genes, n_samples = expression.shape
    expression = np.nan_to_num(expression, nan=0.0, posinf=0.0, neginf=0.0)

    raise_if_cancelled(cancel_check)
    expr_sub, genes_sub, _ = select_top_variable_genes(expression, gene_names, top_n)
    top_n = len(genes_sub)

    row_dendro_data = None
    row_order = np.arange(top_n)
    if cluster_rows and top_n > 1:
        try:
            raise_if_cancelled(cancel_check)
            if distance == "correlation":
                row_dist = np.nan_to_num(pdist(expr_sub, metric="correlation"), nan=0.0)
            else:
                row_dist = pdist(expr_sub, metric=distance)
            row_link = linkage(row_dist, method=method)
            row_d = dendrogram(row_link, no_plot=True)
            row_order = leaves_list(row_link)
            row_dendro_data = {"icoord": row_d["icoord"], "dcoord": row_d["dcoord"]}
        except Exception as exc:
            logger.warning("Row clustering failed: %s", exc)

    col_dendro_data = None
    col_order = np.arange(n_samples)
    if cluster_cols and n_samples > 1:
        try:
            raise_if_cancelled(cancel_check)
            col_dist = np.nan_to_num(pdist(expr_sub.T, metric=distance), nan=0.0)
            col_link = linkage(col_dist, method=method)
            col_d = dendrogram(col_link, no_plot=True)
            col_order = leaves_list(col_link)
            col_dendro_data = {"icoord": col_d["icoord"], "dcoord": col_d["dcoord"]}
        except Exception as exc:
            logger.warning("Column clustering failed: %s", exc)

    ordered = expr_sub[row_order][:, col_order]
    raise_if_cancelled(cancel_check)

    row_means = np.mean(ordered, axis=1, keepdims=True)
    row_stds = np.std(ordered, axis=1, keepdims=True)
    row_stds = np.maximum(row_stds, 1e-6)
    z_scored = (ordered - row_means) / row_stds
    z_scored = np.nan_to_num(z_scored, nan=0.0, posinf=ZSCORE_CLIP_LIMIT, neginf=-ZSCORE_CLIP_LIMIT)

    clip_mask = (z_scored > ZSCORE_CLIP_LIMIT) | (z_scored < -ZSCORE_CLIP_LIMIT)
    n_clipped = int(np.count_nonzero(clip_mask))
    z_scored = np.clip(z_scored, -ZSCORE_CLIP_LIMIT, ZSCORE_CLIP_LIMIT)

    warnings = []
    if n_clipped:
        warnings.append(
            f"Clipped {n_clipped} z-score value(s) to +/-{int(ZSCORE_CLIP_LIMIT)} for stable visualization."
        )

    return {
        "z": z_scored.tolist(),
        "x": [sample_names[i] for i in col_order],
        "y": [genes_sub[i] for i in row_order],
        "row_dendrogram": row_dendro_data,
        "col_dendrogram": col_dendro_data,
        "n_total_genes": n_genes,
        "n_shown_genes": top_n,
        "warnings": warnings,
        "zscore_clipped": bool(n_clipped),
        "n_clipped_values": n_clipped,
    }
