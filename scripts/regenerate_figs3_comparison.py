#!/usr/bin/env python3
"""Regenerate Supplementary Figure S3 (comparison with published CMS classifiers)
with the corrected-pipeline values.

Only the MLHeatmap main-pipeline bars changed in the revision:
  - RF+Forward accuracy 0.900 -> 0.898 (main_result.json accuracy)
  - RF+Forward macro AUC 0.9739 -> 0.9727 (mean of roc_data AUCs)
  - External val (GSE39582) accuracy 0.829 -> 0.748 (gse39582_validation_summary.json)
  - External val (GSE39582) macro AUC 0.9612 -> 0.917 (gse39582_validation_summary.json)
The best-internal LightGBM figures (0.912 / 0.9805) come from the 20-combination
screen, which was NOT rerun in the revision (stated in Results), and the published
comparator values are literature figures carried over from the original submission.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "generated" / "analysis"

RED = "#E64B35"   # MLHeatmap (this study)
CYAN = "#4DBBD5"  # published methods


def load_current_values():
    with (ANALYSIS_DIR / "main_result.json").open() as fh:
        main = json.load(fh)
    with (ANALYSIS_DIR / "gse39582_validation_summary.json").open() as fh:
        gse = json.load(fh)
    rf_acc = main["accuracy"]
    rf_macro = float(np.mean([c["auc"] for c in main["roc_data"]]))
    ext_acc = gse["accuracy"]
    ext_macro = gse["macro_auc"]
    return rf_acc, rf_macro, ext_acc, ext_macro


def main():
    rf_acc, rf_macro, ext_acc, ext_macro = load_current_values()
    print(f"RF+Forward acc={rf_acc:.4f} macroAUC={rf_macro:.4f}; External acc={ext_acc:.4f} macroAUC={ext_macro:.4f}")

    # (label, value, is_mlheatmap); literature/retained values from the original Fig S3.
    acc = [
        ("MLHeatmap (LightGBM)", 0.912, True),        # best-internal screen (not rerun)
        ("MLHeatmap (RF+Forward)", round(rf_acc, 3), True),
        ("CMSclassifier (RF)", 0.870, False),
        ("CMSclassifier (SSP)", 0.850, False),
        ("CMSFFPE (RF)", 0.833, False),
        ("CMScaller (NTP)", 0.830, False),
        ("MLHeatmap (External val.)", round(ext_acc, 3), True),
    ]
    auc = [
        ("MLHeatmap (LightGBM)", 0.9805, True),        # best-internal screen (not rerun)
        ("MLHeatmap (RF+Forward)", round(rf_macro, 4), True),
        ("MLHeatmap (External val.)", round(ext_macro, 4), True),
        ("imCMS (DL)", 0.845, False),
    ]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    fig.suptitle("Comparison with Published CMS Classifiers", fontsize=12, y=1.02)

    def draw(ax, data, title, decimals):
        data = sorted(data, key=lambda t: t[1], reverse=True)
        labels = [d[0] for d in data]
        vals = [d[1] for d in data]
        colors = [RED if d[2] else CYAN for d in data]
        y = np.arange(len(data))[::-1]
        ax.barh(y, vals, color=colors, height=0.62, edgecolor="none")
        for yi, v, is_ml in zip(y, vals, [d[2] for d in data]):
            ax.text(v + 0.004, yi, f"{v:.{decimals}f}", va="center", ha="left",
                    fontsize=8.5, fontweight="bold" if is_ml else "normal")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlim(0.70, 1.00)
        ax.set_xlabel(title.split(" ")[-1], fontsize=9)
        ax.set_title(title, fontsize=10)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    axA.text(-0.02, 1.10, "A", transform=axA.transAxes, fontsize=16, fontweight="bold", va="top")
    axB.text(-0.02, 1.10, "B", transform=axB.transAxes, fontsize=16, fontweight="bold", va="top")
    draw(axA, acc, "Classification Accuracy", 3)
    draw(axB, auc, "Classification AUC", 4)

    handles = [plt.Rectangle((0, 0), 1, 1, color=RED), plt.Rectangle((0, 0), 1, 1, color=CYAN)]
    fig.legend(handles, ["MLHeatmap (this study)", "Published methods"],
               loc="lower center", ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.06))

    fig.tight_layout()
    out = ANALYSIS_DIR / "FigureS3_new.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
