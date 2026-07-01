#!/usr/bin/env python3
"""Batch-effect sensitivity analysis for the CRC-CMS cohort.

Uses the TCGA Tissue Source Site (TSS) code embedded in each sample barcode
as a real technical/site batch proxy (TCGA aggregates sequencing across many
contributing centers; TSS is the standard site covariate used in TCGA batch
QC). Tests:

1. Confounding: is TSS (batch) associated with CMS group (would bias fold
   splits/results if batches cluster by subtype)?
2. Batch-associated expression variance: how much of the top-PC variance in
   normalized expression is explained by TSS vs. by CMS group?
3. A PCA figure colored by TSS and by CMS group, for the supplement.

Usage: python scripts/run_batch_effect_analysis.py [--cohort-dir DIR] [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from mlheatmap.core.gene_mapping import map_gene_ids  # noqa: E402
from mlheatmap.core.normalization import deseq2_normalize  # noqa: E402


def eta_squared_one_way(values: np.ndarray, groups: np.ndarray) -> float:
    """Proportion of variance in `values` explained by categorical `groups`."""
    grand_mean = np.mean(values)
    ss_total = np.sum((values - grand_mean) ** 2)
    ss_between = 0.0
    for g in np.unique(groups):
        mask = groups == g
        group_mean = np.mean(values[mask])
        ss_between += mask.sum() * (group_mean - grand_mean) ** 2
    return float(ss_between / ss_total) if ss_total > 0 else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-dir", type=Path, default=PROJECT_ROOT / "generated" / "crc_cms_public")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "generated" / "analysis")
    parser.add_argument("--min-tss-n", type=int, default=10, help="Minimum samples per TSS to include in tests")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = args.cohort_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    meta_path = args.cohort_dir / "tcga_crc_cms_gold_metadata.tsv"

    print("Loading cohort + metadata...")
    df = pd.read_csv(counts_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.set_index("sample_id").loc[df.columns]  # align exactly to matrix column order
    meta["tss"] = meta["tcga_barcode"].str.split("-").str[1]

    print("Gene mapping + global DESeq2-like normalization (same as main analysis)...")
    gene_ids = df.index.astype(str).tolist()
    mapping, _unmapped = map_gene_ids(gene_ids, "human", "auto")
    df2 = df.copy()
    df2.index = df2.index.astype(str)
    df2["_symbol"] = df2.index.map(lambda g: mapping.get(str(g), str(g)))
    mapped = df2.groupby("_symbol").sum()
    raw_counts = mapped.to_numpy(dtype=np.float64)
    normalized = deseq2_normalize(raw_counts)  # genes x samples

    cms = meta["CMS_gold"].to_numpy()
    tss = meta["tss"].to_numpy()

    # --- 1. Confounding: TSS vs CMS group ---
    tss_counts = pd.Series(tss).value_counts()
    keep_tss = set(tss_counts[tss_counts >= args.min_tss_n].index)
    mask_keep = np.array([t in keep_tss for t in tss])
    contingency = pd.crosstab(tss[mask_keep], cms[mask_keep])
    chi2, p_chi2, dof, _expected = stats.chi2_contingency(contingency)
    cramers_v = float(np.sqrt(chi2 / (contingency.to_numpy().sum() * (min(contingency.shape) - 1))))

    print(f"\nTSS categories with >= {args.min_tss_n} samples: {sorted(keep_tss)} ({mask_keep.sum()} samples)")
    print(f"Chi-square test (TSS x CMS): chi2={chi2:.2f}, dof={dof}, p={p_chi2:.4g}, Cramer's V={cramers_v:.3f}")

    # --- 2. Variance explained by batch vs. by biology on top PCs ---
    # Use the 2,000 most variable genes (same order of magnitude as the app's
    # internal prefilter) so this mirrors what feeds the classifier.
    var_per_gene = np.var(normalized, axis=1)
    top_idx = np.argsort(var_per_gene)[-2000:]
    X = normalized[top_idx, :].T  # samples x genes
    X_centered = X - X.mean(axis=0, keepdims=True)
    # SVD-based PCA (no sklearn dependency needed for this step)
    U, S, _Vt = np.linalg.svd(X_centered, full_matrices=False)
    n_pcs = 10
    pcs = (U[:, :n_pcs] * S[:n_pcs])
    explained_var_ratio = (S**2 / np.sum(S**2))[:n_pcs]

    pc_variance_table = []
    for pc_i in range(n_pcs):
        eta2_cms = eta_squared_one_way(pcs[:, pc_i], cms)
        eta2_tss_all = eta_squared_one_way(pcs[mask_keep, pc_i], tss[mask_keep])
        pc_variance_table.append(
            {
                "pc": pc_i + 1,
                "explained_var_ratio": float(explained_var_ratio[pc_i]),
                "eta2_by_cms_group": eta2_cms,
                "eta2_by_tss_batch": eta2_tss_all,
            }
        )
        print(
            f"  PC{pc_i + 1}: var_ratio={explained_var_ratio[pc_i]:.3f}  "
            f"eta2(CMS)={eta2_cms:.3f}  eta2(TSS)={eta2_tss_all:.3f}"
        )

    # --- 3. Figure: PCA colored by TSS and by CMS ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    cms_palette = {"CMS1": "#E64B35", "CMS2": "#4DBBD5", "CMS3": "#00A087", "CMS4": "#3C5488"}
    for label, color in cms_palette.items():
        m = cms == label
        axes[0].scatter(pcs[m, 0], pcs[m, 1], s=14, alpha=0.7, color=color, label=label)
    axes[0].set_xlabel(f"PC1 ({explained_var_ratio[0] * 100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({explained_var_ratio[1] * 100:.1f}%)")
    axes[0].set_title("Colored by CMS group")
    axes[0].legend(frameon=False, fontsize=8)

    top_tss = list(tss_counts.head(8).index)
    cmap = plt.get_cmap("tab10")
    for i, t in enumerate(top_tss):
        m = tss == t
        axes[1].scatter(pcs[m, 0], pcs[m, 1], s=14, alpha=0.7, color=cmap(i % 10), label=f"{t} (n={m.sum()})")
    m_other = ~np.isin(tss, top_tss)
    axes[1].scatter(pcs[m_other, 0], pcs[m_other, 1], s=10, alpha=0.35, color="lightgray", label="other TSS")
    axes[1].set_xlabel(f"PC1 ({explained_var_ratio[0] * 100:.1f}%)")
    axes[1].set_ylabel(f"PC2 ({explained_var_ratio[1] * 100:.1f}%)")
    axes[1].set_title("Colored by TCGA Tissue Source Site (batch proxy)")
    axes[1].legend(frameon=False, fontsize=6.5, ncol=2)

    fig.suptitle("Batch-effect check: TCGA CRC-CMS cohort, top 2,000 variable genes")
    fig.tight_layout()
    fig_path = args.out_dir / "batch_effect_pca.png"
    fig.savefig(fig_path, dpi=200)
    print(f"\nWrote {fig_path}")

    summary = {
        "n_samples": int(len(cms)),
        "n_distinct_tss": int(pd.Series(tss).nunique()),
        "min_tss_n_for_test": args.min_tss_n,
        "tss_categories_tested": sorted(keep_tss),
        "n_samples_in_confounding_test": int(mask_keep.sum()),
        "confounding_chi2": float(chi2),
        "confounding_dof": int(dof),
        "confounding_p_value": float(p_chi2),
        "confounding_cramers_v": cramers_v,
        "pc_variance_table": pc_variance_table,
        "figure": str(fig_path),
    }
    summary_path = args.out_dir / "batch_effect_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
