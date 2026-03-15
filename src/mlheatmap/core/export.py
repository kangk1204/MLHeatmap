"""Server-side export helpers."""

from __future__ import annotations

import io
from typing import Any

from mlheatmap import __version__


def build_results_metadata(session) -> dict[str, Any]:
    """Build a reproducibility snapshot for the current session."""
    return {
        "app": {
            "name": "MLHeatmap",
            "version": __version__,
        },
        "session": {
            "id": session.id,
            "species": session.species,
            "id_type": session.id_type,
            "norm_method": session.norm_method,
            "deg_effect_size_basis": session.deg_effect_size_basis,
            "n_genes": len(session.gene_names),
            "n_samples": len(session.sample_names),
            "excluded_samples": list(session.excluded_samples),
        },
        "groups": session.groups,
        "metadata": session.metadata,
    }


def _flatten_metadata(data: Any, prefix: str = "") -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if isinstance(data, dict):
        for key, value in data.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(_flatten_metadata(value, child_prefix))
        return rows
    if isinstance(data, list):
        for idx, value in enumerate(data):
            child_prefix = f"{prefix}[{idx}]"
            rows.extend(_flatten_metadata(value, child_prefix))
        return rows
    rows.append({"field": prefix, "value": data})
    return rows


def export_results_excel(session) -> bytes:
    """Export all current analysis results and metadata to an Excel workbook."""
    import pandas as pd

    output = io.BytesIO()
    metadata = build_results_metadata(session)

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        pd.DataFrame(_flatten_metadata(metadata)).to_excel(writer, sheet_name="Metadata", index=False)

        if session.normalized is not None:
            df = pd.DataFrame(
                session.normalized,
                index=session.gene_names,
                columns=session.sample_names,
            )
            df.index.name = "Gene"
            df.to_excel(writer, sheet_name="Normalized Expression")

        if session.biomarker_results:
            genes_df = pd.DataFrame(session.biomarker_results["top_genes"])
            genes_df.to_excel(writer, sheet_name="Biomarker Genes", index=False)
            optimal_combo = session.biomarker_results.get("optimal_combo")
            if optimal_combo:
                pd.DataFrame(optimal_combo.get("auc_curve", [])).to_excel(
                    writer,
                    sheet_name="Panel AUC Curve",
                    index=False,
                )
                pd.DataFrame(optimal_combo.get("selection_frequency", [])).to_excel(
                    writer,
                    sheet_name="Panel Consensus",
                    index=False,
                )

        if session.deg_results:
            deg_data = session.deg_results["results"]
            deg_df = pd.DataFrame(deg_data)
            col_order = [
                "gene",
                "log2fc",
                "pvalue",
                "adj_pvalue",
                "direction",
                "mean_g1",
                "mean_g2",
                "neg_log10_p",
                "gene_idx",
            ]
            col_order = [column for column in col_order if column in deg_df.columns]
            deg_df = deg_df[col_order]
            comparison_group = session.deg_results.get("comparison_group", "")
            reference_group = session.deg_results.get("reference_group", "")
            if comparison_group and reference_group:
                deg_df = deg_df.rename(
                    columns={
                        "mean_g1": f"mean_{comparison_group}",
                        "mean_g2": f"mean_{reference_group}",
                    }
                )
            deg_df.to_excel(writer, sheet_name="DEG Results", index=False)

        if session.groups:
            rows = []
            for group, samples in session.groups.items():
                for sample in samples:
                    rows.append({"Sample": sample, "Group": group})
            pd.DataFrame(rows).to_excel(writer, sheet_name="Groups", index=False)

    return output.getvalue()
