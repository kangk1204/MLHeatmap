#!/usr/bin/env python3
"""Second independent external validation: CMCBSN Korean CRC RNA-seq cohort.

Transfers the TCGA-derived compact gene panel to an independent Korean CRC
RNA-seq cohort (Lee et al. 2024, Mol Cells 47(3):100033; expression + CMS
labels from Zenodo doi:10.5281/zenodo.8333650), following the same
z-score-then-classify transfer recipe already used for the GSE39582
microarray cohort in this manuscript, but with both cohorts DESeq2-like
normalized first (both are RNA-seq counts here, unlike the microarray case).

Usage:
  python scripts/run_cmcbsn_validation.py \
      --tcga-cohort-dir generated/crc_cms_public \
      --cmcbsn-dir generated/cmcbsn/sources \
      --panel-genes-from generated/analysis/main_result.json \
      --out-dir generated/analysis
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelBinarizer

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from mlheatmap.core.biomarker import _build_model  # noqa: E402
from mlheatmap.core.gene_mapping import map_gene_ids  # noqa: E402
from mlheatmap.core.normalization import deseq2_normalize  # noqa: E402

EXPECTED_MD5 = {
    "CMCBSN_expectedcount_381.txt": "8d3a19c707648e12f62256e60eaebf40",
    "Supplementary table5.xlsx": "ea5181bdc4823fe87cc5f634342821ef",
}


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_cmcbsn_cms_labels(supp_table5_path: Path) -> pd.Series:
    """Sample-id (dash format, matching the count-matrix header) -> CMS label,
    tumor samples only, dropping the ~14% with no confident CMS call (as
    reported by the original study)."""
    wb = openpyxl.load_workbook(supp_table5_path, data_only=True)
    ws = wb["information"]
    rows = list(ws.iter_rows(values_only=True))
    records = []
    for row in rows[1:]:
        if row[0] is None:
            continue
        dot_id = str(row[0])
        dash_id = dot_id.replace(".", "-")
        cms = row[-1]
        records.append((dash_id, cms))
    frame = pd.DataFrame(records, columns=["sample_id", "CMS"])
    frame = frame[frame["CMS"].isin(["CMS1", "CMS2", "CMS3", "CMS4"])]
    return frame.set_index("sample_id")["CMS"]


def load_and_map_counts(counts_path: Path, sep: str = "\t") -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(counts_path, sep=sep, index_col=0)
    gene_ids = df.index.astype(str).tolist()
    mapping, unmapped = map_gene_ids(gene_ids, "human", "auto")
    df = df.copy()
    df.index = df.index.astype(str)
    df["_symbol"] = df.index.map(lambda g: mapping.get(str(g), str(g)))
    mapped = df.groupby("_symbol").sum()
    mapped.index.name = None
    meta = {"mapped_count": len(mapping), "unmapped_count": len(unmapped), "total": len(gene_ids)}
    return mapped, meta


def zscore_rows(mat: np.ndarray) -> np.ndarray:
    mean = mat.mean(axis=1, keepdims=True)
    std = mat.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (mat - mean) / std


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tcga-cohort-dir", type=Path, default=PROJECT_ROOT / "generated" / "crc_cms_public")
    parser.add_argument("--cmcbsn-dir", type=Path, default=PROJECT_ROOT / "generated" / "cmcbsn" / "sources")
    parser.add_argument("--panel-genes-from", type=Path, default=PROJECT_ROOT / "generated" / "analysis" / "main_result.json")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "generated" / "analysis")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = args.cmcbsn_dir / "CMCBSN_expectedcount_381.txt"
    labels_path = args.cmcbsn_dir / "supp_table5.xlsx"
    for path, name in [
        (counts_path, "CMCBSN_expectedcount_381.txt"),
        (labels_path, "Supplementary table5.xlsx"),
    ]:
        expected = EXPECTED_MD5.get(name)
        if expected is not None:
            actual = md5_file(path)
            status = "OK" if actual == expected else "MISMATCH"
            print(f"md5 {name}: {actual} ({status} vs pinned {expected})")

    with args.panel_genes_from.open() as fh:
        main_result = json.load(fh)
    panel_genes = main_result["optimal_combo"]["best_genes"]
    print(f"TCGA-derived panel ({len(panel_genes)} genes): {panel_genes}")

    print("\nLoading + gene-mapping CMCBSN cohort...")
    cmcbsn_counts, cmcbsn_map_meta = load_and_map_counts(counts_path)
    print(f"  {cmcbsn_map_meta}, {cmcbsn_counts.shape[0]} genes x {cmcbsn_counts.shape[1]} samples (all, N+T)")

    cms_labels = load_cmcbsn_cms_labels(labels_path)
    print(f"  {len(cms_labels)} tumor samples with a confident CMS call: {cms_labels.value_counts().to_dict()}")

    missing_samples = [s for s in cms_labels.index if s not in cmcbsn_counts.columns]
    if missing_samples:
        print(f"  WARNING: {len(missing_samples)} labeled samples not found in count matrix: {missing_samples[:5]}...")
    cms_labels = cms_labels[[s for s in cms_labels.index if s in cmcbsn_counts.columns]]

    cmcbsn_tumor = cmcbsn_counts[cms_labels.index]
    print(f"  Final CMCBSN validation set: {cmcbsn_tumor.shape[1]} tumor samples")

    panel_in_cmcbsn = [g for g in panel_genes if g in cmcbsn_tumor.index]
    panel_missing_cmcbsn = [g for g in panel_genes if g not in cmcbsn_tumor.index]
    print(f"\n  Panel genes represented in CMCBSN: {len(panel_in_cmcbsn)}/{len(panel_genes)} (missing: {panel_missing_cmcbsn})")

    print("\nLoading + gene-mapping TCGA training cohort...")
    tcga_counts_path = args.tcga_cohort_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    tcga_counts, tcga_map_meta = load_and_map_counts(tcga_counts_path)
    tcga_meta = pd.read_csv(args.tcga_cohort_dir / "tcga_crc_cms_gold_metadata.tsv", sep="\t").set_index("sample_id")
    tcga_meta = tcga_meta.loc[tcga_counts.columns]
    print(f"  {tcga_counts.shape[0]} genes x {tcga_counts.shape[1]} samples")

    panel_in_tcga = [g for g in panel_in_cmcbsn if g in tcga_counts.index]
    if len(panel_in_tcga) != len(panel_in_cmcbsn):
        print(f"  WARNING: {set(panel_in_cmcbsn) - set(panel_in_tcga)} in CMCBSN but not in TCGA (unexpected)")
    used_genes = panel_in_tcga
    print(f"\nUsing {len(used_genes)}/{len(panel_genes)} panel genes present in both cohorts.")

    print("\nNormalizing both cohorts (DESeq2-like VST, each cohort independently)...")
    tcga_norm_full = deseq2_normalize(tcga_counts.to_numpy(dtype=np.float64))
    tcga_norm = pd.DataFrame(tcga_norm_full, index=tcga_counts.index, columns=tcga_counts.columns)
    cmcbsn_norm_full = deseq2_normalize(cmcbsn_tumor.to_numpy(dtype=np.float64))
    cmcbsn_norm = pd.DataFrame(cmcbsn_norm_full, index=cmcbsn_tumor.index, columns=cmcbsn_tumor.columns)

    X_train_raw = tcga_norm.loc[used_genes].to_numpy(dtype=np.float64)  # genes x samples
    X_test_raw = cmcbsn_norm.loc[used_genes].to_numpy(dtype=np.float64)

    X_train = zscore_rows(X_train_raw).T  # samples x genes
    X_test = zscore_rows(X_test_raw).T

    y_train = tcga_meta["CMS_gold"].to_numpy()
    y_test = cms_labels.to_numpy()

    print("\nTraining Random Forest (same _build_model config as the main pipeline, n_estimators=500, random_state=42) on full TCGA panel-gene matrix...")
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
    print(f"\nCMCBSN external validation: accuracy={accuracy:.4f}, macro_auc={macro_auc:.4f}")
    print(confusion)

    summary = {
        "cohort": "CMCBSN (Lee et al. 2024, Mol Cells; Zenodo 10.5281/zenodo.8333650)",
        "n_tumor_samples_total_with_cms_call": int(len(cms_labels)),
        "n_samples_used": int(len(y_test)),
        "cms_distribution": cms_labels.value_counts().to_dict(),
        "panel_genes_total": len(panel_genes),
        "panel_genes_used": used_genes,
        "panel_genes_missing_in_cmcbsn": panel_missing_cmcbsn,
        "n_panel_genes_used": len(used_genes),
        "accuracy": accuracy,
        "macro_auc": macro_auc,
        "confusion_matrix": confusion.to_dict(),
        "gene_mapping_cmcbsn": cmcbsn_map_meta,
        "gene_mapping_tcga": tcga_map_meta,
    }
    out_path = args.out_dir / "cmcbsn_validation_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
