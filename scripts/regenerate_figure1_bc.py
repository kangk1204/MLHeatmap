#!/usr/bin/env python3
"""Regenerate Figure 1 panels B (per-class AUC) and C (SHAP heatmap) from the
corrected-pipeline main_result.json, then composite with panels A and D
cropped from the original submitted Figure1.png (UI screenshot and the
20-combination screen, neither of which changed) into a new Figure1.png.

Also regenerates Figure S2 (compact-panel gene bars) to reflect the new panel.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "generated" / "analysis"
MANUSCRIPT_DIR = Path("C:/Projects/09_MLHeatmap_Revision/101_Manuscript")

CMS_COLORS = {"CMS1": "#E64B35", "CMS2": "#4DBBD5", "CMS3": "#00A087", "CMS4": "#3C5488"}
CMS_DISPLAY = {"CMS1": "MSI-Immune", "CMS2": "Canonical", "CMS3": "Metabolic", "CMS4": "Mesenchymal"}


def make_panel_b(result: dict, out_path: Path):
    fig, ax = plt.subplots(figsize=(3.6, 4.6))
    fig.text(0.01, 0.98, "B", fontsize=28, fontweight="bold", ha="left", va="top")
    groups = [c["group"] for c in result["roc_data"]]
    aucs = [c["auc"] for c in result["roc_data"]]
    stds = [c["std"] for c in result["roc_data"]]
    colors = [CMS_COLORS[g] for g in groups]
    x = np.arange(len(groups))
    ax.bar(x, aucs, yerr=stds, color=colors, capsize=4, width=0.65, edgecolor="none")
    for xi, v in zip(x, aucs):
        ax.text(xi, v + 0.012, f"{v:.3f}", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([CMS_DISPLAY[g] for g in groups], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("AUC (out-of-fold)", fontsize=10)
    ax.set_ylim(0.5, 1.05)
    ax.set_title("Classification performance", fontsize=11)
    macro_auc = float(np.mean(aucs))
    ax.text(
        0.5, 0.15,
        f"Overall accuracy: {result['accuracy'] * 100:.1f}%\nMacro AUC: {macro_auc:.3f}",
        transform=ax.transAxes, ha="center", va="center", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"),
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, transparent=True)
    plt.close(fig)
    return macro_auc


def make_panel_c(result: dict, out_path: Path):
    genes = [g["gene"] for g in result["shap_plot_data"]]
    shap_abs = np.array([np.mean(g["values"]) for g in result["shap_plot_data"]])
    expr = np.array([g["expression"] for g in result["shap_plot_data"]])  # genes x samples
    labels = np.array(result["sample_labels"])

    order = np.argsort(labels, kind="stable")
    labels_sorted = labels[order]
    expr_sorted = expr[:, order]

    z = (expr_sorted - expr_sorted.mean(axis=1, keepdims=True)) / (expr_sorted.std(axis=1, keepdims=True) + 1e-9)
    z = np.clip(z, -2.5, 2.5)

    fig = plt.figure(figsize=(6.2, 4.8))
    fig.text(0.01, 0.98, "C", fontsize=28, fontweight="bold", ha="left", va="top")
    gs = gridspec.GridSpec(2, 2, width_ratios=[10, 2.2], height_ratios=[0.5, 10], hspace=0.03, wspace=0.15)

    ax_anno = fig.add_subplot(gs[0, 0])
    for cms in ["CMS1", "CMS2", "CMS3", "CMS4"]:
        idx = np.where(labels_sorted == cms)[0]
        if len(idx) == 0:
            continue
        ax_anno.axvspan(idx.min(), idx.max() + 1, color=CMS_COLORS[cms])
        ax_anno.text((idx.min() + idx.max() + 1) / 2, 0.5, cms, ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax_anno.set_xlim(0, z.shape[1])
    ax_anno.set_ylim(0, 1)
    ax_anno.axis("off")

    ax_heat = fig.add_subplot(gs[1, 0])
    im = ax_heat.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="none")
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    ax_heat.set_ylabel(f"Top {len(genes)} genes", fontsize=10)

    ax_bar = fig.add_subplot(gs[1, 1])
    ax_bar.barh(range(len(genes)), shap_abs[::-1], color="#6A51A3", height=0.7)
    ax_bar.set_ylim(-0.5, len(genes) - 0.5)
    ax_bar.invert_yaxis()
    ax_bar.set_yticks([])
    ax_bar.set_title("SHAP\nvalue", fontsize=9)
    for spine in ("top", "right", "left"):
        ax_bar.spines[spine].set_visible(False)

    fig.savefig(out_path, dpi=300, transparent=True)
    plt.close(fig)


def make_fig_s2(result: dict, out_path: Path):
    genes = result["optimal_combo"]["best_genes"]
    n = len(genes)
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, n - 1)) for i in range(n)]
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    y = np.arange(n)[::-1]
    ax.barh(y, [1] * n, color=colors, height=0.75)
    for yi, g, rank in zip(y, genes, range(1, n + 1)):
        ax.text(1.02, yi, f"#{rank}", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlim(0, 1.15)
    ax.set_xticks([])
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.set_title(f"Compact Panel: {n} genes (AUC = {result['optimal_combo']['best_auc']:.3f})", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, transparent=False, facecolor="white")
    plt.close(fig)


def crop_panels_a_d(out_a: Path, out_d: Path):
    """Crop panels A (UI screenshot) and D (20-combination screen) from the
    original submitted Figure1.png (2208x675); these panels did not change in
    the revision. Panel D is cropped from x=1470 so its "D" label (originally
    spanning ~x=1500-1560) is fully included with margin, rather than sliced."""
    # Crop source is the ORIGINAL submitted Figure1.png (old B/C panels, but the
    # A UI screenshot and the D 20-combination screen are reused unchanged). A
    # preserved copy lives beside the analysis outputs so this stays correct even
    # after the submission-folder Figure1.png is replaced with the new composite.
    preserved = ANALYSIS_DIR / "Figure1_source_original.png"
    source = preserved if preserved.exists() else (MANUSCRIPT_DIR / "Figure1.png")
    original = Image.open(source).convert("RGB")
    w, h = original.size  # 2208 x 675
    original.crop((0, 0, 408, h)).save(out_a)
    original.crop((1470, 0, w, h)).save(out_d)
    print(f"Cropped panel A [0:408] -> {out_a.name}, panel D [1470:{w}] -> {out_d.name}")


def flatten_to_white(im: Image.Image) -> Image.Image:
    """Composite an RGBA image onto an opaque white background (avoids
    alpha-blending artifacts, e.g. thin text disappearing, on later resize)."""
    if im.mode != "RGBA":
        return im.convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    bg.paste(im, (0, 0), im)
    return bg


def composite_figure1(panel_b_path: Path, panel_c_path: Path, out_path: Path):
    panelA = flatten_to_white(Image.open(ANALYSIS_DIR / "fig1_panelA.png"))
    panelD = flatten_to_white(Image.open(ANALYSIS_DIR / "fig1_panelD.png"))
    panelB = flatten_to_white(Image.open(panel_b_path))
    panelC = flatten_to_white(Image.open(panel_c_path))

    h = 675
    def resize_to_h(im, target_h):
        w = int(im.width * target_h / im.height)
        return im.resize((w, target_h), Image.LANCZOS)

    panelA_r = resize_to_h(panelA, h)
    panelD_r = resize_to_h(panelD, h)
    panelB_r = resize_to_h(panelB, h)
    panelC_r = resize_to_h(panelC, h)

    total_w = panelA_r.width + panelB_r.width + panelC_r.width + panelD_r.width
    canvas = Image.new("RGB", (total_w, h), (255, 255, 255))
    x = 0
    labels = ["A", "B", "C", "D"]
    for panel, label in zip((panelA_r, panelB_r, panelC_r, panelD_r), labels):
        canvas.paste(panel, (x, 0))
        print(f"panel {label}: x={x}, w={panel.width}")
        x += panel.width
    canvas.save(out_path)
    print(f"Wrote composited {out_path} ({canvas.size})")


def main():
    with (ANALYSIS_DIR / "main_result.json").open() as fh:
        result = json.load(fh)

    crop_panels_a_d(ANALYSIS_DIR / "fig1_panelA.png", ANALYSIS_DIR / "fig1_panelD.png")

    panel_b_path = ANALYSIS_DIR / "fig1_panelB_new.png"
    panel_c_path = ANALYSIS_DIR / "fig1_panelC_new.png"
    macro_auc = make_panel_b(result, panel_b_path)
    print(f"Panel B written, macro_auc={macro_auc:.4f}")
    make_panel_c(result, panel_c_path)
    print("Panel C written")

    composite_figure1(panel_b_path, panel_c_path, ANALYSIS_DIR / "Figure1_new.png")

    make_fig_s2(result, ANALYSIS_DIR / "FigureS2_new.png")
    print("Fig S2 written")


if __name__ == "__main__":
    main()
