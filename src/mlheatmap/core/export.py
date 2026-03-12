"""Export functions for images and Excel."""

import io
import numpy as np


def _check_kaleido():
    """Check that kaleido is installed for image export."""
    try:
        import kaleido  # noqa: F401
    except ImportError:
        raise ImportError(
            "kaleido is required for image export. "
            "Install it with: pip install kaleido"
        )


def export_heatmap_image(session, fmt="png", dpi=300) -> bytes:
    """Export heatmap as PNG or SVG using plotly."""
    import plotly.graph_objects as go
    _check_kaleido()

    data = session.heatmap_data
    if not data:
        raise ValueError("No heatmap data. Generate heatmap first.")

    cs = getattr(session, "heatmap_color_scale", "RdBu_r") or "RdBu_r"
    fig = go.Figure(data=go.Heatmap(
        z=data["z"],
        x=data["x"],
        y=data["y"],
        colorscale=cs,
        zmid=0,
        colorbar=dict(title="Z-score"),
    ))

    fig.update_layout(
        width=max(800, len(data["x"]) * 25 + 200),
        height=max(600, len(data["y"]) * 3 + 200),
        font=dict(family="Arial", size=10),
        xaxis=dict(tickangle=-45),
        margin=dict(l=120, r=60, t=40, b=120),
    )

    if fmt == "svg":
        data = fig.to_image(format="svg", scale=2)
        return data.encode() if isinstance(data, str) else data
    return fig.to_image(format="png", scale=dpi / 100)


def export_shap_image(session, fmt="png", dpi=300) -> bytes:
    """Export SHAP importance bar chart."""
    import plotly.graph_objects as go
    _check_kaleido()

    results = session.biomarker_results
    if not results:
        raise ValueError("No biomarker results. Run analysis first.")

    genes = [g["symbol"] for g in results["top_genes"]]
    values = [g["shap_mean_abs"] for g in results["top_genes"]]

    fig = go.Figure(data=go.Bar(
        x=list(reversed(values)),
        y=list(reversed(genes)),
        orientation="h",
        marker=dict(
            color=list(reversed(values)),
            colorscale="Viridis",
        ),
    ))

    fig.update_layout(
        title="Top Biomarker Genes (SHAP Importance)",
        xaxis_title="Mean |SHAP value|",
        width=800,
        height=max(400, len(genes) * 25 + 100),
        font=dict(family="Arial", size=12),
        margin=dict(l=120, r=40, t=60, b=60),
    )

    if fmt == "svg":
        return fig.to_image(format="svg", scale=2)
    return fig.to_image(format="png", scale=dpi / 100)


def export_auc_image(session, fmt="png", dpi=300) -> bytes:
    """Export ROC/AUC curve."""
    import plotly.graph_objects as go
    _check_kaleido()

    results = session.biomarker_results
    if not results or "roc_data" not in results:
        raise ValueError("No AUC data. Run analysis first.")

    colors = ["#3b82f6", "#ef4444", "#10b981", "#f59e0b", "#8b5cf6"]
    fig = go.Figure()

    for i, curve in enumerate(results["roc_data"]):
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(
            x=curve["fpr"],
            y=curve["tpr"],
            name=f'{curve["group"]} (AUC={curve["auc"]:.3f}±{curve["std"]:.3f})',
            line=dict(color=color, width=2),
        ))

    fig.add_trace(go.Scatter(
        x=[0, 1], y=[0, 1],
        line=dict(color="gray", dash="dash", width=1),
        showlegend=False,
    ))

    fig.update_layout(
        title="ROC Curve (Out-of-Fold)",
        xaxis_title="False Positive Rate",
        yaxis_title="True Positive Rate",
        width=700,
        height=600,
        font=dict(family="Arial", size=12),
        legend=dict(x=0.5, y=0.05),
    )

    if fmt == "svg":
        return fig.to_image(format="svg", scale=2)
    return fig.to_image(format="png", scale=dpi / 100)


def export_volcano_image(session, fmt="png", dpi=300) -> bytes:
    """Export Volcano plot as publication-quality figure using matplotlib."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    deg = session.deg_results
    if not deg or "results" not in deg:
        raise ValueError("No DEG results. Run DEG analysis first.")

    results = deg["results"]
    thresholds = deg["thresholds"]
    fc_thresh = thresholds["log2fc"]
    p_thresh = thresholds["pvalue"]
    neg_log10_thresh = -np.log10(p_thresh) if p_thresh > 0 else 2
    pvalue_type = deg.get("pvalue_type", "fdr")
    p_label = "P-value" if pvalue_type == "raw" else "FDR"

    # Separate by direction
    up = [r for r in results if r["direction"] == "up"]
    down = [r for r in results if r["direction"] == "down"]
    ns = [r for r in results if r["direction"] == "ns"]

    # Dark theme matching the app
    BG = "#06060b"
    BG_CARD = "#0f1117"
    TEXT = "#e2e8f0"
    TEXT_SEC = "#94a3b8"
    GRID = "#1e293b"

    fig, ax = plt.subplots(figsize=(9, 7), facecolor=BG)
    ax.set_facecolor(BG_CARD)

    # NS points (small gray)
    if ns:
        ax.scatter(
            [r["log2fc"] for r in ns],
            [r["neg_log10_p"] for r in ns],
            s=6, c="#4b5563", alpha=0.35, edgecolors="none",
            label=f"NS ({len(ns):,})", rasterized=True, zorder=1,
        )
    # Down points (blue)
    if down:
        ax.scatter(
            [r["log2fc"] for r in down],
            [r["neg_log10_p"] for r in down],
            s=30, c="#3b82f6", alpha=0.8, edgecolors="none",
            label=f"Down ({len(down)})", zorder=2,
        )
    # Up points (red)
    if up:
        ax.scatter(
            [r["log2fc"] for r in up],
            [r["neg_log10_p"] for r in up],
            s=30, c="#ef4444", alpha=0.8, edgecolors="none",
            label=f"Up ({len(up)})", zorder=2,
        )

    # Top gene labels (up to 5 per direction, spread vertically to reduce overlap)
    top_up = sorted(
        up,
        key=lambda r: r["pvalue"] if pvalue_type == "raw" else r.get("adj_pvalue", r["pvalue"]),
    )[:5]
    top_down = sorted(
        down,
        key=lambda r: r["pvalue"] if pvalue_type == "raw" else r.get("adj_pvalue", r["pvalue"]),
    )[:5]
    label_genes = top_up + top_down
    offsets_y = [8, 20, 32, 44, 56, 8, 20, 32, 44, 56]  # stagger vertically
    for i, r in enumerate(label_genes):
        ax.annotate(
            r["gene"],
            (r["log2fc"], r["neg_log10_p"]),
            fontsize=7.5, color=TEXT, alpha=0.9, fontweight="500",
            textcoords="offset points", xytext=(0, offsets_y[i % len(offsets_y)]),
            ha="center", va="bottom",
            arrowprops=dict(arrowstyle="-", color=TEXT_SEC, alpha=0.4, lw=0.5),
        )

    # Threshold lines
    all_fc = [r["log2fc"] for r in results]
    max_fc = max(abs(min(all_fc)), abs(max(all_fc)), fc_thresh + 0.5) if all_fc else 3
    ax.axvline(fc_thresh, color="#475569", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(-fc_thresh, color="#475569", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(neg_log10_thresh, color="#475569", linestyle="--", linewidth=0.8, alpha=0.6)

    # Styling
    ax.set_xlabel("log₂ Fold Change", fontsize=12, color=TEXT_SEC, labelpad=8)
    ax.set_ylabel(f"-log₁₀({p_label})", fontsize=12, color=TEXT_SEC, labelpad=8)
    ax.set_xlim(-max_fc * 1.15, max_fc * 1.15)
    ax.tick_params(colors=TEXT_SEC, labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.spines["left"].set_color(GRID)
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.15, color=GRID)

    legend = ax.legend(
        loc="upper left", fontsize=9, frameon=True,
        facecolor=BG_CARD, edgecolor=GRID, labelcolor=TEXT,
    )
    legend.get_frame().set_alpha(0.9)

    n_total = len(results)
    n_sig = len(up) + len(down)
    method = deg.get("method", "wilcoxon")
    method_label = "Wilcoxon" if method == "wilcoxon" else "Welch's t-test"
    fig.suptitle(
        f"Volcano Plot — {n_sig:,} DEGs / {n_total:,} genes ({method_label})",
        fontsize=14, color=TEXT, y=1.02, fontweight="500",
    )
    fig.text(
        0.99, 0.01,
        f"|log₂FC| > {fc_thresh}  &  {p_label} < {p_thresh}",
        fontsize=8, color=TEXT_SEC, alpha=0.7, ha="right",
        transform=fig.transFigure,
    )

    buf = io.BytesIO()
    fig.savefig(
        buf, format=fmt, dpi=dpi, bbox_inches="tight",
        facecolor=BG, edgecolor="none", pad_inches=0.2,
    )
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def export_results_excel(session) -> bytes:
    """Export all results to Excel."""
    import pandas as pd

    output = io.BytesIO()

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        # Sheet 1: Normalized expression
        if session.normalized is not None:
            df = pd.DataFrame(
                session.normalized,
                index=session.gene_names,
                columns=session.sample_names,
            )
            df.index.name = "Gene"
            df.to_excel(writer, sheet_name="Normalized Expression")

        # Sheet 2: Biomarker results
        if session.biomarker_results:
            genes_df = pd.DataFrame(session.biomarker_results["top_genes"])
            genes_df.to_excel(writer, sheet_name="Biomarker Genes", index=False)

        # Sheet 3: DEG results
        if hasattr(session, "deg_results") and session.deg_results:
            deg_data = session.deg_results["results"]
            deg_df = pd.DataFrame(deg_data)
            # Reorder columns for clarity
            col_order = ["gene", "log2fc", "pvalue", "adj_pvalue", "direction",
                         "mean_g1", "mean_g2", "neg_log10_p", "gene_idx"]
            col_order = [c for c in col_order if c in deg_df.columns]
            deg_df = deg_df[col_order]
            deg_df.to_excel(writer, sheet_name="DEG Results", index=False)

        # Sheet 4: Group assignments
        if session.groups:
            rows = []
            for group, samples in session.groups.items():
                for s in samples:
                    rows.append({"Sample": s, "Group": group})
            pd.DataFrame(rows).to_excel(writer, sheet_name="Groups", index=False)

    return output.getvalue()
