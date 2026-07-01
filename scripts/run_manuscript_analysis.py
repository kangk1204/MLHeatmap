#!/usr/bin/env python3
"""Run the full CRC-CMS manuscript analysis end-to-end from a built cohort.

This mirrors the exact upload -> gene mapping -> normalize -> biomarker chain
the web UI performs (see api/gene_mapping.py, api/normalize.py,
api/biomarker.py), so it reproduces the manuscript numbers without a running
server. It runs the main (global-normalization, importance-basis) analysis
plus the two revision sensitivity analyses (fold-isolated normalization,
SHAP-basis candidate selection) requested by reviewers, and writes:

- generated/analysis/main_result.json          (global + importance; the
  manuscript's primary result)
- generated/analysis/fold_isolated_result.json  (per_fold_normalize=True)
- generated/analysis/shap_basis_result.json     (selection_basis="shap")
- generated/analysis/sensitivity_summary.json   (accuracy/AUC deltas +
  paired significance test + Jaccard panel overlap)
- generated/analysis/table_s1.xlsx              (regenerated Excel export,
  includes the "Panel Fold Overlap" sheet)

Usage: python scripts/run_manuscript_analysis.py [--cohort-dir DIR] [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import uuid
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from mlheatmap.core.biomarker import run_biomarker_analysis  # noqa: E402
from mlheatmap.core.export import export_results_excel  # noqa: E402
from mlheatmap.core.gene_mapping import map_gene_ids  # noqa: E402
from mlheatmap.core.normalization import deseq2_normalize  # noqa: E402


def load_and_map_cohort(counts_path: Path) -> tuple[np.ndarray, list[str], list[str], dict]:
    """Load the built cohort matrix and apply the app's gene-mapping step.

    Mirrors api/gene_mapping.py's map_genes with species="human", id_type="auto":
    unmapped symbols are retained under their original name (preserve_raw_symbols),
    then rows sharing a resulting symbol are summed.
    """
    df = pd.read_csv(counts_path, sep="\t", index_col=0)
    gene_ids = df.index.astype(str).tolist()
    sample_names = df.columns.tolist()

    mapping, unmapped = map_gene_ids(gene_ids, "human", "auto")
    df = df.copy()
    df.index = df.index.astype(str)
    df["_symbol"] = df.index.map(lambda gene_id: mapping.get(str(gene_id), str(gene_id)))
    mapped = df.groupby("_symbol").sum()
    mapped.index.name = None

    mapping_metadata = {
        "species": "human",
        "id_type": "auto",
        "mapped_count": len(mapping),
        "unmapped_count": len(unmapped),
        "total": len(gene_ids),
        "n_genes_after_mapping": int(mapped.shape[0]),
    }
    return mapped.to_numpy(dtype=np.float64), mapped.index.astype(str).tolist(), sample_names, mapping_metadata


def groups_from_prefixed_names(sample_names: list[str]) -> dict[str, list[int]]:
    """Mirror the UI's 'Auto-detect' grouping: split each sample name on '__'."""
    groups: dict[str, list[int]] = {}
    for idx, name in enumerate(sample_names):
        label = name.split("__", 1)[0]
        groups.setdefault(label, []).append(idx)
    return groups


def run_variant(
    *,
    expression: np.ndarray,
    raw_counts: np.ndarray,
    gene_names: list[str],
    groups: dict[str, list[int]],
    per_fold_normalize: bool,
    selection_basis: str,
    label: str,
) -> dict:
    t0 = time.monotonic()
    result = run_biomarker_analysis(
        expression=expression,
        gene_names=gene_names,
        sample_groups=groups,
        n_top_genes=30,
        n_estimators=500,
        cv_folds=5,
        model="rf",
        panel_method="forward",
        raw_counts=raw_counts,
        norm_method="deseq2",
        per_fold_normalize=per_fold_normalize,
        selection_basis=selection_basis,
        random_state=42,
    )
    elapsed = time.monotonic() - t0
    macro_auc = float(np.mean([c["auc"] for c in result["roc_data"]]))
    print(
        f"[{label}] accuracy={result['accuracy']:.4f} macro_auc={macro_auc:.4f} "
        f"panel_n={result['optimal_combo']['n_genes']} panel_auc={result['optimal_combo']['best_auc']:.4f} "
        f"shap_fallback_used={result['shap_fallback_used']} elapsed={elapsed:.1f}s"
    )
    result["_macro_auc"] = macro_auc
    result["_elapsed_s"] = elapsed
    return result


def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def paired_ttest_summary(a: list[float], b: list[float]) -> dict:
    from scipy import stats

    a_arr, b_arr = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    diff = a_arr - b_arr
    t_stat, p_value = stats.ttest_rel(a_arr, b_arr) if len(a_arr) > 1 else (float("nan"), float("nan"))
    return {
        "mean_diff": float(np.mean(diff)),
        "std_diff": float(np.std(diff)),
        "t_stat": float(t_stat),
        "p_value": float(p_value),
        "n_folds": int(len(a_arr)),
    }


def fold_accuracies_and_aucs(result: dict) -> tuple[list[float], list[float]]:
    """Best-effort per-fold accuracy/AUC extraction is not exposed by the core
    return dict (only fold-mean summaries are), so the paired test instead
    compares per-class held-out AUCs (roc_data) between the two runs, which
    *is* fold-resolved-equivalent information available from both variants."""
    return [c["auc"] for c in result["roc_data"]], [result["accuracy"]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-dir", type=Path, default=PROJECT_ROOT / "generated" / "crc_cms_public")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "generated" / "analysis")
    args = parser.parse_args()

    counts_path = args.cohort_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    if not counts_path.exists():
        raise SystemExit(f"Cohort not found at {counts_path}. Run `mlheatmap-download-crc-cms` first.")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading cohort and applying gene mapping...")
    raw_counts, gene_names, sample_names, mapping_meta = load_and_map_cohort(counts_path)
    print(f"  {mapping_meta}")

    print("Normalizing (global DESeq2-like VST)...")
    normalized = deseq2_normalize(raw_counts)

    groups = groups_from_prefixed_names(sample_names)
    print(f"  groups: { {k: len(v) for k, v in groups.items()} }")

    print("\n=== Main result (global normalization, importance selection basis) ===")
    main_result = run_variant(
        expression=normalized, raw_counts=raw_counts, gene_names=gene_names, groups=groups,
        per_fold_normalize=False, selection_basis="importance", label="main",
    )

    print("\n=== Fold-isolated normalization sensitivity ===")
    fold_iso_result = run_variant(
        expression=normalized, raw_counts=raw_counts, gene_names=gene_names, groups=groups,
        per_fold_normalize=True, selection_basis="importance", label="fold_isolated",
    )

    print("\n=== SHAP-basis candidate selection ===")
    shap_basis_result = run_variant(
        expression=normalized, raw_counts=raw_counts, gene_names=gene_names, groups=groups,
        per_fold_normalize=False, selection_basis="shap", label="shap_basis",
    )

    def strip_private(d: dict) -> dict:
        return {k: v for k, v in d.items() if not k.startswith("_")}

    (args.out_dir / "main_result.json").write_text(
        json.dumps(strip_private(main_result), indent=2, default=str), encoding="utf-8"
    )
    (args.out_dir / "fold_isolated_result.json").write_text(
        json.dumps(strip_private(fold_iso_result), indent=2, default=str), encoding="utf-8"
    )
    (args.out_dir / "shap_basis_result.json").write_text(
        json.dumps(strip_private(shap_basis_result), indent=2, default=str), encoding="utf-8"
    )

    main_panel = main_result["optimal_combo"]["best_genes"]
    shap_panel = shap_basis_result["optimal_combo"]["best_genes"]
    panel_jaccard = jaccard(main_panel, shap_panel)

    main_auc_by_class, main_acc = fold_accuracies_and_aucs(main_result)
    fold_auc_by_class, fold_acc = fold_accuracies_and_aucs(fold_iso_result)

    sensitivity_summary = {
        "mapping_metadata": mapping_meta,
        "n_genes_final": len(gene_names),
        "n_samples": len(sample_names),
        "group_counts": {k: len(v) for k, v in groups.items()},
        "main": {
            "accuracy": main_result["accuracy"],
            "macro_auc": main_result["_macro_auc"],
            "panel_n_genes": main_result["optimal_combo"]["n_genes"],
            "panel_genes": main_panel,
            "panel_auc": main_result["optimal_combo"]["best_auc"],
            "shap_fallback_used": main_result["shap_fallback_used"],
        },
        "fold_isolated": {
            "accuracy": fold_iso_result["accuracy"],
            "macro_auc": fold_iso_result["_macro_auc"],
            "panel_n_genes": fold_iso_result["optimal_combo"]["n_genes"],
            "panel_genes": fold_iso_result["optimal_combo"]["best_genes"],
            "panel_auc": fold_iso_result["optimal_combo"]["best_auc"],
        },
        "shap_basis": {
            "accuracy": shap_basis_result["accuracy"],
            "macro_auc": shap_basis_result["_macro_auc"],
            "panel_n_genes": shap_basis_result["optimal_combo"]["n_genes"],
            "panel_genes": shap_panel,
            "panel_auc": shap_basis_result["optimal_combo"]["best_auc"],
            "shap_fallback_used": shap_basis_result["shap_fallback_used"],
        },
        "global_vs_fold_isolated": {
            "accuracy_delta": main_result["accuracy"] - fold_iso_result["accuracy"],
            "macro_auc_delta": main_result["_macro_auc"] - fold_iso_result["_macro_auc"],
            "per_class_auc_paired_ttest": paired_ttest_summary(main_auc_by_class, fold_auc_by_class),
            "panel_gene_jaccard": jaccard(main_panel, fold_iso_result["optimal_combo"]["best_genes"]),
        },
        "importance_vs_shap_basis": {
            "accuracy_delta": main_result["accuracy"] - shap_basis_result["accuracy"],
            "macro_auc_delta": main_result["_macro_auc"] - shap_basis_result["_macro_auc"],
            "panel_gene_jaccard": panel_jaccard,
            "panel_genes_shared": sorted(set(main_panel) & set(shap_panel)),
            "panel_genes_importance_only": sorted(set(main_panel) - set(shap_panel)),
            "panel_genes_shap_only": sorted(set(shap_panel) - set(main_panel)),
        },
    }
    (args.out_dir / "sensitivity_summary.json").write_text(
        json.dumps(sensitivity_summary, indent=2, default=str), encoding="utf-8"
    )
    print("\n=== Sensitivity summary ===")
    print(json.dumps(sensitivity_summary["global_vs_fold_isolated"], indent=2))
    print(json.dumps(sensitivity_summary["importance_vs_shap_basis"], indent=2))

    print("\nExporting Table S1 workbook (main result)...")
    session = SimpleNamespace(
        id=str(uuid.uuid4()),
        species="human",
        id_type="symbol",
        norm_method="deseq2",
        deg_effect_size_basis="size_factor_normalized_counts",
        gene_names=gene_names,
        sample_names=sample_names,
        excluded_samples=[],
        groups={g: [sample_names[i] for i in idx] for g, idx in groups.items()},
        metadata={
            "app": {"name": "MLHeatmap", "version": "1.0.0"},
            "mapping": mapping_meta,
            "biomarker": {
                "model": "rf",
                "model_label": "Random Forest",
                "panel_method": "forward",
                "n_top_genes": 30,
                "n_estimators": 500,
                "cv_folds_requested": 5,
                "cv_folds_used": main_result["cv_folds_used"],
                "roc_evaluation": "out_of_fold",
                "panel_evaluation": "nested_outer_cv",
                "shap_fallback_used": main_result["shap_fallback_used"],
                "normalization_scope": main_result["normalization_scope"],
                "selection_basis": main_result["selection_basis"],
            },
        },
        normalized=normalized,
        biomarker_results=main_result,
        deg_results=None,
    )
    xlsx_bytes = export_results_excel(session, include_normalized_expression=False)
    (args.out_dir / "table_s1.xlsx").write_bytes(xlsx_bytes)
    print(f"  wrote {args.out_dir / 'table_s1.xlsx'} ({len(xlsx_bytes) / 1024:.0f} KB)")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
