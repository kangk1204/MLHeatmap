"""Server-side export helpers."""

import io


def export_results_excel(session) -> bytes:
    """Export all current analysis results to an Excel workbook."""
    import pandas as pd

    output = io.BytesIO()

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
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

        if getattr(session, "deg_results", None):
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
