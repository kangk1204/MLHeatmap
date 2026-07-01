#!/usr/bin/env python3
"""UpSet-style fold-overlap figure: which panel genes were selected in which
of the 5 outer CV folds, sorted by selection frequency."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "generated" / "analysis"


def main():
    with (ANALYSIS_DIR / "main_result.json").open() as fh:
        result = json.load(fh)
    oc = result["optimal_combo"]
    fold_panels = oc["fold_panels"]
    sel_freq = oc["selection_frequency"]

    genes_sorted = [d["gene"] for d in sel_freq]
    freqs = {d["gene"]: d["frequency"] for d in sel_freq}
    n_folds = len(fold_panels)
    n_genes = len(genes_sorted)

    membership = np.zeros((n_genes, n_folds), dtype=bool)
    for j, fp in enumerate(fold_panels):
        gene_set = set(fp["genes"])
        for i, g in enumerate(genes_sorted):
            membership[i, j] = g in gene_set

    fig, (ax_mat, ax_bar) = plt.subplots(
        1, 2, figsize=(7.5, 0.28 * n_genes + 1.2), gridspec_kw={"width_ratios": [3, 1.3]}, sharey=True
    )

    for i in range(n_genes):
        for j in range(n_folds):
            color = "#3C5488" if membership[i, j] else "#E8E8E8"
            ax_mat.scatter(j, i, s=90, color=color, zorder=3)
        row_present = np.where(membership[i])[0]
        if len(row_present) > 1:
            ax_mat.plot([row_present.min(), row_present.max()], [i, i], color="#3C5488", lw=1.5, zorder=2, alpha=0.5)

    ax_mat.set_xticks(range(n_folds))
    ax_mat.set_xticklabels([f"Fold {j + 1}" for j in range(n_folds)], fontsize=9)
    ax_mat.set_yticks(range(n_genes))
    ax_mat.set_yticklabels(genes_sorted, fontsize=8)
    ax_mat.invert_yaxis()
    ax_mat.set_xlim(-0.5, n_folds - 0.5)
    ax_mat.set_ylim(n_genes - 0.5, -0.5)
    for spine in ("top", "right", "left", "bottom"):
        ax_mat.spines[spine].set_visible(False)
    ax_mat.set_title("Fold membership", fontsize=10)

    ax_bar.barh(range(n_genes), [freqs[g] * n_folds for g in genes_sorted], color="#4DBBD5", height=0.6)
    ax_bar.invert_yaxis()
    ax_bar.set_xlim(0, n_folds + 0.3)
    ax_bar.set_xticks(range(0, n_folds + 1))
    ax_bar.set_xlabel("Folds selected in\n(of 5)", fontsize=9)
    for spine in ("top", "right"):
        ax_bar.spines[spine].set_visible(False)

    fig.suptitle("Consensus compact-panel gene selection across outer CV folds", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = ANALYSIS_DIR / "fig_fold_overlap_upset.png"
    fig.savefig(out_path, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_path}")

    n_all_1 = sum(1 for d in sel_freq if d["frequency"] == 1.0)
    print(f"Genes selected in all {n_folds} folds: {n_all_1} ({[d['gene'] for d in sel_freq if d['frequency'] == 1.0]})")
    print(f"Total distinct genes appearing across any fold: {n_genes}")


if __name__ == "__main__":
    main()
