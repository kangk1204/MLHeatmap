#!/usr/bin/env python3
"""External validation: transfer the TCGA-derived compact panel to GSE39582
(Affymetrix GPL570, Marisa et al. 2013), the same recipe described in the
manuscript (z-score each cohort independently at the gene level, transfer a
classifier trained on the full TCGA panel-gene matrix), but using CMScaller
-derived CMS labels since the original CRCSC network-consensus label source
for this cohort is no longer reachable (see Methods note).

Usage: python scripts/run_gse39582_validation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelBinarizer

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from mlheatmap.core.biomarker import _build_model  # noqa: E402
from mlheatmap.core.gene_mapping import map_gene_ids  # noqa: E402
from mlheatmap.core.normalization import deseq2_normalize  # noqa: E402


def zscore_rows(mat: np.ndarray) -> np.ndarray:
    mean = mat.mean(axis=1, keepdims=True)
    std = mat.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (mat - mean) / std


def main() -> int:
    tcga_dir = PROJECT_ROOT / "generated" / "crc_cms_public"
    gse_dir = PROJECT_ROOT / "generated" / "gse39582"
    out_dir = PROJECT_ROOT / "generated" / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    with (out_dir / "main_result.json").open() as fh:
        main_result = json.load(fh)
    panel_genes = main_result["optimal_combo"]["best_genes"]
    print(f"TCGA-derived panel ({len(panel_genes)} genes): {panel_genes}")

    print("\nLoading GSE39582 tumor expression + CMScaller labels...")
    gse_expr = pd.read_csv(gse_dir / "gse39582_tumor_gene_expr.tsv.gz", sep="\t", index_col=0)
    labels = pd.read_csv(gse_dir / "gse39582_cmscaller_labels.csv").set_index("sample")
    labels = labels[labels["cms_label"].isin(["CMS1", "CMS2", "CMS3", "CMS4"])]
    print(f"  {gse_expr.shape[0]} genes x {gse_expr.shape[1]} samples; {len(labels)} with a confident CMS call")

    gse_tumor = gse_expr[labels.index]

    panel_in_gse = [g for g in panel_genes if g in gse_tumor.index]
    panel_missing_gse = [g for g in panel_genes if g not in gse_tumor.index]
    print(f"\n  Panel genes represented on GPL570: {len(panel_in_gse)}/{len(panel_genes)} (missing: {panel_missing_gse})")

    print("\nLoading + gene-mapping TCGA training cohort...")
    tcga_counts_path = tcga_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    tcga_counts = pd.read_csv(tcga_counts_path, sep="\t", index_col=0)
    gene_ids = tcga_counts.index.astype(str).tolist()
    mapping, _unmapped = map_gene_ids(gene_ids, "human", "auto")
    tcga_counts.index = tcga_counts.index.astype(str)
    tcga_counts["_symbol"] = tcga_counts.index.map(lambda g: mapping.get(str(g), str(g)))
    tcga_counts = tcga_counts.groupby("_symbol").sum()
    tcga_counts.index.name = None
    tcga_meta = pd.read_csv(tcga_dir / "tcga_crc_cms_gold_metadata.tsv", sep="\t").set_index("sample_id")
    tcga_meta = tcga_meta.loc[tcga_counts.columns]
    print(f"  {tcga_counts.shape[0]} genes x {tcga_counts.shape[1]} samples")

    used_genes = [g for g in panel_in_gse if g in tcga_counts.index]
    print(f"\nUsing {len(used_genes)}/{len(panel_genes)} panel genes present in both cohorts.")

    print("\nNormalizing TCGA (DESeq2-like VST); GSE39582 already log2 RMA-normalized (microarray)...")
    tcga_norm_full = deseq2_normalize(tcga_counts.to_numpy(dtype=np.float64))
    tcga_norm = pd.DataFrame(tcga_norm_full, index=tcga_counts.index, columns=tcga_counts.columns)

    X_train_raw = tcga_norm.loc[used_genes].to_numpy(dtype=np.float64)
    X_test_raw = gse_tumor.loc[used_genes].to_numpy(dtype=np.float64)

    X_train = zscore_rows(X_train_raw).T
    X_test = zscore_rows(X_test_raw).T

    y_train = tcga_meta["CMS_gold"].to_numpy()
    y_test = labels["cms_label"].to_numpy()

    print("\nTraining Random Forest (same _build_model config as main pipeline) on full TCGA panel-gene matrix...")
    clf, _use_shap_tree, _needs_scaling = _build_model("rf", n_estimators=500, n_samples=len(y_train), random_state=42)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    accuracy = float(np.mean(y_pred == y_test))

    lb = LabelBinarizer()
    lb.fit(sorted(set(y_train)))
    y_test_bin = lb.transform(y_test)
    y_proba = clf.predict_proba(X_test)
    class_order = list(clf.classes_)
    y_proba_aligned = np.column_stack([y_proba[:, class_order.index(c)] for c in lb.classes_])
    macro_auc = float(roc_auc_score(y_test_bin, y_proba_aligned, average="macro", multi_class="ovr"))

    confusion = pd.crosstab(pd.Series(y_test, name="gold"), pd.Series(y_pred, name="predicted"))
    print(f"\nGSE39582 external validation: accuracy={accuracy:.4f}, macro_auc={macro_auc:.4f}")
    print(confusion)

    summary = {
        "cohort": "GSE39582 (Marisa et al. 2013, Affymetrix GPL570)",
        "label_source": "CMScaller-derived (original CRCSC network-consensus source no longer reachable)",
        "n_tumor_samples_total": int(gse_expr.shape[1]),
        "n_samples_with_confident_cms": int(len(labels)),
        "cms_distribution": labels["cms_label"].value_counts().to_dict(),
        "panel_genes_total": len(panel_genes),
        "panel_genes_used": used_genes,
        "panel_genes_missing_on_gpl570": panel_missing_gse,
        "n_panel_genes_used": len(used_genes),
        "accuracy": accuracy,
        "macro_auc": macro_auc,
        "confusion_matrix": confusion.to_dict(),
    }
    out_path = out_dir / "gse39582_validation_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
