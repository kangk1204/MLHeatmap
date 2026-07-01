#!/usr/bin/env python3
"""Build the real head-to-head comparison table: CRCSC gold label vs. CMScaller
vs. MLHeatmap, sample-by-sample, on the identical 511-sample TCGA cohort.

Requires:
  - generated/analysis/main_result.json (from run_manuscript_analysis.py; must
    include sample_indices + oof_predicted_labels, added to run_biomarker_analysis
    in this revision)
  - generated/analysis/cmscaller_predictions.csv (from run_cmscaller_benchmark.R)
  - generated/crc_cms_public/tcga_crc_cms_gold_counts.tsv.gz (for sample order)

Usage: python scripts/compare_cmscaller_vs_mlheatmap.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    analysis_dir = PROJECT_ROOT / "generated" / "analysis"
    cohort_dir = PROJECT_ROOT / "generated" / "crc_cms_public"

    with (analysis_dir / "main_result.json").open() as fh:
        main_result = json.load(fh)

    counts_path = cohort_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    sample_names = pd.read_csv(counts_path, sep="\t", index_col=0, nrows=0).columns.tolist()

    sample_indices = main_result["sample_indices"]
    true_labels = main_result["sample_labels"]
    mlheatmap_pred = main_result["oof_predicted_labels"]

    mlheatmap_df = pd.DataFrame(
        {
            "sample": [sample_names[i] for i in sample_indices],
            "gold_label": true_labels,
            "mlheatmap_prediction": mlheatmap_pred,
        }
    ).set_index("sample")

    cmscaller_df = pd.read_csv(analysis_dir / "cmscaller_predictions.csv").set_index("sample")

    merged = mlheatmap_df.join(cmscaller_df[["cmscaller_prediction"]], how="inner")
    assert len(merged) == len(mlheatmap_df), "sample set mismatch between MLHeatmap and CMScaller runs"

    merged["mlheatmap_correct"] = merged["gold_label"] == merged["mlheatmap_prediction"]
    merged["cmscaller_confident"] = merged["cmscaller_prediction"].notna()
    merged["cmscaller_correct"] = merged["cmscaller_confident"] & (
        merged["gold_label"] == merged["cmscaller_prediction"]
    )
    merged["both_confident_agree"] = merged["cmscaller_confident"] & (
        merged["mlheatmap_prediction"] == merged["cmscaller_prediction"]
    )

    n = len(merged)
    mlheatmap_acc = merged["mlheatmap_correct"].mean()
    cmscaller_confident = merged[merged["cmscaller_confident"]]
    cmscaller_acc_confident = cmscaller_confident["cmscaller_correct"].mean()
    agreement_confident = merged.loc[merged["cmscaller_confident"], "both_confident_agree"].mean()

    print(f"n={n} samples (identical cohort, identical fold-blind CRCSC gold labels)")
    print(f"MLHeatmap (out-of-fold) accuracy vs. gold: {mlheatmap_acc:.4f}")
    print(
        f"CMScaller accuracy vs. gold (confident calls only, "
        f"{merged['cmscaller_confident'].sum()}/{n}): {cmscaller_acc_confident:.4f}"
    )
    print(f"MLHeatmap <-> CMScaller agreement (confident calls only): {agreement_confident:.4f}")

    confusion_ml = pd.crosstab(merged["gold_label"], merged["mlheatmap_prediction"])
    confusion_cms_ml = pd.crosstab(
        cmscaller_confident["cmscaller_prediction"], cmscaller_confident["mlheatmap_prediction"]
    )
    print("\nMLHeatmap confusion (gold rows x predicted cols):")
    print(confusion_ml)
    print("\nCMScaller vs MLHeatmap agreement matrix (CMScaller rows x MLHeatmap cols, confident CMScaller calls):")
    print(confusion_cms_ml)

    out = {
        "n_samples": int(n),
        "mlheatmap_accuracy": float(mlheatmap_acc),
        "cmscaller_n_confident": int(merged["cmscaller_confident"].sum()),
        "cmscaller_accuracy_confident": float(cmscaller_acc_confident),
        "mlheatmap_cmscaller_agreement_confident": float(agreement_confident),
        "mlheatmap_confusion": confusion_ml.to_dict(),
        "mlheatmap_vs_cmscaller_confusion": confusion_cms_ml.to_dict(),
    }
    out_path = analysis_dir / "head_to_head_summary.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    merged.to_csv(analysis_dir / "head_to_head_per_sample.csv")
    print(f"\nWrote {out_path} and head_to_head_per_sample.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
