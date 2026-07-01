#!/usr/bin/env python3
"""Rerun just the main (global normalization, importance selection basis) variant
and overwrite generated/analysis/main_result.json — used to pick up code changes
made to run_biomarker_analysis after the original 3-variant run had already
started (e.g. the sample_indices/oof_predicted_labels fields added mid-run)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_manuscript_analysis import load_and_map_cohort, groups_from_prefixed_names, run_variant  # noqa: E402
from mlheatmap.core.normalization import deseq2_normalize  # noqa: E402

cohort_dir = PROJECT_ROOT / "generated" / "crc_cms_public"
out_dir = PROJECT_ROOT / "generated" / "analysis"
counts_path = cohort_dir / "tcga_crc_cms_gold_counts.tsv.gz"

print("Loading cohort and applying gene mapping...")
raw_counts, gene_names, sample_names, mapping_meta = load_and_map_cohort(counts_path)
normalized = deseq2_normalize(raw_counts)
groups = groups_from_prefixed_names(sample_names)

print("\n=== Main result (rerun to pick up sample_indices/oof_predicted_labels) ===")
main_result = run_variant(
    expression=normalized, raw_counts=raw_counts, gene_names=gene_names, groups=groups,
    per_fold_normalize=False, selection_basis="importance", label="main",
)

def strip_private(d):
    return {k: v for k, v in d.items() if not k.startswith("_")}

out_path = out_dir / "main_result.json"
out_path.write_text(json.dumps(strip_private(main_result), indent=2, default=str), encoding="utf-8")
print(f"Overwrote {out_path}")
print("Has sample_indices:", "sample_indices" in main_result)
print("Has oof_predicted_labels:", "oof_predicted_labels" in main_result)
