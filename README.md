# MLHeatmap

Interactive RNA-seq heatmap and biomarker discovery tool.

## Install

```bash
pip install mlheatmap
```

## Run

```bash
mlheatmap
```

Auto-opens browser at http://127.0.0.1:8765

## Features

- Upload bulk RNA-seq raw count matrices (CSV/TSV/Excel)
- Auto gene ID mapping (Human/Mouse) to gene symbols
- Normalization: DESeq2-like VST, CPM+Log2, Log2(count+1)
- Interactive clustered heatmap with dendrograms
- Drag-and-drop sample group assignment
- Biomarker discovery via Random Forest + SHAP
- Cross-validated ROC/AUC curves
- High-resolution image (PNG/SVG) and Excel export
