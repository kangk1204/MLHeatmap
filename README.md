# MLHeatmap

Interactive RNA-seq heatmap and biomarker discovery tool for bulk count matrices.

## Quickstart

### Native install

From a local clone of this repository:

```bash
pip install .
mlheatmap
```

The app starts on `http://127.0.0.1:8765`.

To run the test suite from a local clone:

```bash
pip install .[test]
pytest
```

### One-step installers

Windows:

```bat
install-windows.cmd
```

Ubuntu:

```bash
bash ./install-ubuntu.sh
```

macOS Apple Silicon:

```bash
bash ./install-macos.sh
```

Each installer:

- Reuses Python `3.11` or `3.12`, or bootstraps one when needed
- Creates or repairs a local `.venv`
- Installs MLHeatmap as a regular package
- Runs `mlheatmap --self-check`
- Starts the app locally

After installation, restart with:

```bash
bash ./run-ubuntu.sh --no-browser
```

Windows equivalent:

```bat
run-windows.cmd --no-browser
```

### Docker

Docker is the recommended path when you need the optional `xgboost` and `lightgbm` models:

```bash
docker build -t mlheatmap .
docker run --rm -p 8765:8765 mlheatmap
```

## Input Rules

MLHeatmap expects a matrix with:

- Gene IDs in the first column
- Sample names in the header row
- Numeric count values in every remaining cell

Accepted formats:

- `.csv`
- `.tsv`
- `.txt`
- `.xlsx`
- `.xls`

Fail-closed validation:

- Fully empty rows or columns are dropped
- Any remaining missing value, text value, mixed-type cell, `NaN`, or `inf` causes upload failure
- The upload API returns `invalid_cell_count`, `invalid_examples`, and `invalid_columns` for malformed files

Automatic preprocessing:

- All-zero genes are removed
- Low-expression genes are filtered with the default rule: count `>= 10` in at least `max(2, 20% of samples)` samples

## Methodology And Limitations

### Normalization

- `DESeq2-like VST`: median-of-ratios normalization with a VST-style transform
- `TPM / CPM + log2`: linear TPM when gene lengths are available, otherwise CPM fallback, then `log2(x + 1)`
- `log2(count + 1)`: direct log transform of counts

`DESeq2-like VST` is a normalization label only. It is not the DESeq2 differential-expression test.

### Biomarker analysis

- SHAP ranking and ROC curves are computed with outer-fold held-out evaluation
- Compact panel selection (`forward`, `lasso`, `stability`, `mrmr`) now uses nested outer CV
- `optimal_combo.best_auc` is the mean held-out AUC from the outer folds
- `optimal_combo.selection_frequency` reports how often each consensus panel gene was selected across folds

Paper-safe interpretation:

- Use `roc_data` and nested panel AUCs for manuscript performance claims
- Do not reuse results produced before this hardening pass

### DEG analysis

DEG is exploratory in this app.

- Statistical tests: Wilcoxon rank-sum or Welch's t-test
- Testing is performed on normalized expression
- `log2fc` is always computed from a linear-scale abundance basis:
  - raw counts for `log2`
  - TPM/CPM abundance for `tpm`
  - size-factor-normalized counts for `deseq2`

Paper-safe interpretation:

- Describe this as exploratory DEG on normalized data
- Do not describe these results as "DESeq2 differential expression" unless you run an external DESeq2 workflow

## Paper Reproducibility

Use the non-interactive reproduction script after any code or method change:

```bash
python scripts/paper_reproduce.py \
  --input tests/data/human_ensembl_12samples.csv \
  --groups-json path/to/groups.json \
  --output-dir reproduction_out \
  --species human \
  --id-type auto \
  --normalize deseq2 \
  --model rf \
  --panel-method forward
```

`groups.json` must be a JSON object:

```json
{
  "Control": ["Sample1", "Sample2"],
  "Case": ["Sample3", "Sample4"]
}
```

The script writes:

- `biomarker.json`
- `deg.json`
- `metadata.json`
- `results.xlsx`

`results.xlsx` includes a `Metadata` sheet with provenance for:

- upload filtering
- gene mapping
- normalization
- DEG settings
- biomarker settings
- current package version

## Manuscript Checklist

Before submission:

- Re-run all manuscript figures and tables with `scripts/paper_reproduce.py`
- Replace any legacy volcano table or DEG result created before the effect-size fix
- Replace any legacy optimal panel AUC created before nested outer-CV evaluation
- Describe DEG as exploratory Wilcoxon/Welch testing on normalized data
- Do not write "DESeq2 analysis" when referring only to the app's normalization mode
- Record the repository commit and Python version used for the final rerun

## Supported Platforms

- Windows 11
- Ubuntu
- macOS Apple Silicon
- Python 3.11-3.12
