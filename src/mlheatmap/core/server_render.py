"""Server-side heatmap rendering using matplotlib/seaborn."""

import io
import logging
import math

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist

from mlheatmap.core.clustering import select_top_variable_genes

logger = logging.getLogger(__name__)

# App color palette (matches heatmap.js GROUP_COLORS)
GROUP_COLORS = ["#3b82f6", "#ef4444", "#10b981", "#f59e0b", "#8b5cf6",
                "#ec4899", "#14b8a6", "#f97316", "#6366f1", "#84cc16"]

# Dark theme colors (matches app CSS variables)
BG_PRIMARY = "#06060b"
BG_CARD = "#0f1117"
TEXT_PRIMARY = "#e2e8f0"
TEXT_SECONDARY = "#94a3b8"
DENDRO_COLOR = "#475569"

# Plotly → matplotlib colorscale mapping
CMAP_MAP = {
    "RdBu_r": "RdBu_r",
    "RdYlBu_r": "RdYlBu_r",
    "Viridis": "viridis",
    "Inferno": "inferno",
    "Plasma": "plasma",
}

# Row clustering memory threshold (~1GB pdist limit)
MAX_CLUSTER_GENES = 15000


def render_heatmap_image(
    expression: np.ndarray,
    gene_names: list[str],
    sample_names: list[str],
    groups: dict | None = None,
    top_n: int = 500,
    distance: str = "correlation",
    method: str = "average",
    color_scale: str = "RdBu_r",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    fmt: str = "png",
    dpi: int = 150,
) -> bytes:
    """Render a publication-quality clustered heatmap as an image.

    Uses seaborn.clustermap with dark theme styling matching the app.
    Handles large gene sets (up to 60,000) with adaptive clustering.
    """
    # Replace non-finite values (from log2(0) = -inf, etc.)
    expression = np.nan_to_num(expression, nan=0.0, posinf=0.0, neginf=0.0)

    # Select top variable genes
    expr_sub, genes_sub, _ = select_top_variable_genes(expression, gene_names, top_n)
    n_genes = len(genes_sub)
    n_samples = len(sample_names)

    logger.info("Server render: %d genes x %d samples", n_genes, n_samples)

    # Z-score rows
    row_means = np.mean(expr_sub, axis=1, keepdims=True)
    row_stds = np.std(expr_sub, axis=1, keepdims=True)
    row_stds = np.maximum(row_stds, 1e-6)
    z_scored = (expr_sub - row_means) / row_stds
    z_scored = np.nan_to_num(z_scored, nan=0.0, posinf=3.0, neginf=-3.0)
    z_scored = np.clip(z_scored, -3, 3)

    # Build DataFrame
    df = pd.DataFrame(z_scored, index=genes_sub, columns=sample_names)

    # Adaptive row clustering
    do_row_cluster = cluster_rows and n_genes > 1
    if do_row_cluster and n_genes > MAX_CLUSTER_GENES:
        logger.info("Skipping row clustering for %d genes (> %d threshold)", n_genes, MAX_CLUSTER_GENES)
        do_row_cluster = False
        # Sort by variance instead (most variable at top)
        variances = np.var(z_scored, axis=1)
        var_order = np.argsort(variances)[::-1]
        df = df.iloc[var_order]

    do_col_cluster = cluster_cols and n_samples > 1

    # Group color bar
    col_colors = None
    group_color_map = {}
    if groups:
        sample_to_group = {}
        for i, (gname, gsamples) in enumerate(groups.items()):
            color = GROUP_COLORS[i % len(GROUP_COLORS)]
            group_color_map[gname] = color
            for s in gsamples:
                sample_to_group[s] = color
        col_color_list = [sample_to_group.get(s, "#333333") for s in sample_names]
        col_colors = pd.Series(col_color_list, index=sample_names, name="Group")

    # Colorscale
    cmap = CMAP_MAP.get(color_scale, "RdBu_r")

    # Distance metric
    dist_metric = distance if distance != "correlation" else "correlation"

    # Figure size
    fig_width = max(10, n_samples * 0.6 + 4)
    fig_height = max(8, min(100, n_genes * 0.015 + 4))

    # Gene label settings
    if n_genes <= 2000:
        yticklabels = True
        gene_fontsize = max(1, min(8, 2000 / n_genes * 3))
    elif n_genes <= 10000:
        step = math.ceil(n_genes / 200)
        yticklabels = step
        gene_fontsize = max(1, min(5, 200 / (n_genes / step)))
    else:
        yticklabels = False
        gene_fontsize = 1

    # Sample labels
    sample_fontsize = max(5, min(11, 500 / n_samples))

    # Dendrogram ratio (use tiny non-zero to avoid seaborn empty bboxes bug)
    dendro_row_ratio = 0.08 if do_row_cluster else 0.001
    dendro_col_ratio = 0.06 if do_col_cluster else 0.001

    try:
        # Set dark theme
        with plt.rc_context({
            "figure.facecolor": BG_PRIMARY,
            "axes.facecolor": BG_CARD,
            "text.color": TEXT_PRIMARY,
            "axes.labelcolor": TEXT_SECONDARY,
            "xtick.color": TEXT_SECONDARY,
            "ytick.color": TEXT_SECONDARY,
        }):
            g = sns.clustermap(
                df,
                method=method,
                metric=dist_metric,
                row_cluster=do_row_cluster,
                col_cluster=do_col_cluster,
                col_colors=col_colors,
                cmap=cmap,
                vmin=-3,
                vmax=3,
                center=0,
                figsize=(fig_width, fig_height),
                dendrogram_ratio=(dendro_col_ratio, dendro_row_ratio),
                cbar_kws={"label": "Z-score", "shrink": 0.5},
                xticklabels=True,
                yticklabels=yticklabels,
                linewidths=0,
                rasterized=True,  # Faster rendering for large matrices
            )

            # Style the figure
            fig = g.fig
            fig.patch.set_facecolor(BG_PRIMARY)

            # Style heatmap axes
            ax_heatmap = g.ax_heatmap
            ax_heatmap.set_facecolor(BG_CARD)
            ax_heatmap.tick_params(axis="x", labelsize=sample_fontsize, rotation=45)
            ax_heatmap.tick_params(axis="y", labelsize=gene_fontsize)

            # Style dendrograms
            if do_row_cluster and hasattr(g, "ax_row_dendrogram"):
                ax_rd = g.ax_row_dendrogram
                ax_rd.set_facecolor(BG_PRIMARY)
                for coll in ax_rd.collections:
                    coll.set_color(DENDRO_COLOR)

            if do_col_cluster and hasattr(g, "ax_col_dendrogram"):
                ax_cd = g.ax_col_dendrogram
                ax_cd.set_facecolor(BG_PRIMARY)
                for coll in ax_cd.collections:
                    coll.set_color(DENDRO_COLOR)

            # Style color bar
            if hasattr(g, "ax_col_colors") and g.ax_col_colors is not None:
                g.ax_col_colors.set_facecolor(BG_PRIMARY)

            # Style colorbar
            cbar = g.cax
            cbar.set_facecolor(BG_PRIMARY)
            cbar.tick_params(colors=TEXT_SECONDARY, labelsize=8)
            cbar.yaxis.label.set_color(TEXT_SECONDARY)

            # Add group legend
            if groups and group_color_map:
                patches = [
                    mpatches.Patch(color=color, label=gname)
                    for gname, color in group_color_map.items()
                ]
                legend = fig.legend(
                    handles=patches,
                    loc="upper right",
                    bbox_to_anchor=(0.98, 0.98),
                    fontsize=max(7, min(10, 100 / len(groups))),
                    frameon=True,
                    facecolor=BG_CARD,
                    edgecolor="#1e293b",
                    labelcolor=TEXT_PRIMARY,
                    title="Groups",
                    title_fontsize=max(8, min(11, 100 / len(groups))),
                )
                legend.get_title().set_color(TEXT_PRIMARY)

            # Add info text
            info_text = f"Top {n_genes:,} genes (most variable)"
            if n_genes > MAX_CLUSTER_GENES:
                info_text += " | Row clustering disabled (large N)"
            fig.text(
                0.02, 0.01, info_text,
                fontsize=7, color=TEXT_SECONDARY, alpha=0.7,
                transform=fig.transFigure,
            )

            # Title
            fig.suptitle(
                f"Clustered Heatmap — {n_genes:,} genes × {n_samples} samples",
                fontsize=12, color=TEXT_PRIMARY, y=0.99, fontweight="500",
            )

            # Render to bytes
            buf = io.BytesIO()
            fig.savefig(
                buf,
                format=fmt,
                dpi=dpi,
                bbox_inches="tight",
                facecolor=BG_PRIMARY,
                edgecolor="none",
                pad_inches=0.3,
            )
            plt.close(fig)
            buf.seek(0)
            return buf.getvalue()

    except Exception as e:
        plt.close("all")
        logger.error("Server render failed: %s", e)
        raise
