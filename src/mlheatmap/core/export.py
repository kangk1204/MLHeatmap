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

    fig = go.Figure(data=go.Heatmap(
        z=data["z"],
        x=data["x"],
        y=data["y"],
        colorscale="RdBu_r",
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
        title="ROC Curve (Cross-Validated)",
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

        # Sheet 3: Group assignments
        if session.groups:
            rows = []
            for group, samples in session.groups.items():
                for s in samples:
                    rows.append({"Sample": s, "Group": group})
            pd.DataFrame(rows).to_excel(writer, sheet_name="Groups", index=False)

    return output.getvalue()
