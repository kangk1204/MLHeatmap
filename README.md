# MLHeatmap

Interactive RNA-seq heatmap and biomarker discovery tool.

## Supported Platforms

- Windows 11
- Ubuntu
- macOS Apple Silicon
- Python 3.11-3.12

## Native Install (Core)

Native install is the recommended path when you want the most stable cross-platform setup.

```bash
pip install mlheatmap
mlheatmap
```

The app starts on `http://127.0.0.1:8765`.

### Windows One-Step Install

From `cmd.exe`, run:

```bat
install-windows.cmd
```

The script reuses an existing Python `3.12` or `3.11`, or bootstraps Python `3.12` automatically with `winget` when no compatible interpreter is present. It creates a local `.venv`, recreates that `.venv` if it was built with the wrong version or became corrupted, installs the app as a regular package, runs a self-check, and starts MLHeatmap. It does not install Python packages into the user's global environment.

After installation, you can start the app again with:

```bat
run-windows.cmd
```

Extra CLI flags still work:

```bat
run-windows.cmd --no-browser
```

### Ubuntu One-Step Install

Run from the repo root:

```bash
bash ./install-ubuntu.sh
```

The script reuses an existing Python `3.12` or `3.11`, or bootstraps one automatically with `apt-get` when no compatible interpreter is present. If the selected interpreter is missing `venv` support, the installer also adds the matching `python3.x-venv` package automatically. It creates a local `.venv`, recreates that `.venv` if it was built with the wrong version or became corrupted, installs the app as a regular package, runs a self-check, and starts MLHeatmap. It does not install Python packages into the system Python environment.

After installation, you can start the app again with:

```bash
bash ./run-ubuntu.sh
```

Extra CLI flags still work:

```bash
bash ./run-ubuntu.sh --no-browser
```

### macOS Apple Silicon One-Step Install

Run from the repo root:

```bash
bash ./install-macos.sh
```

The script reuses an existing Python `3.12` or `3.11`, or bootstraps Python `3.12` automatically with Homebrew when no compatible interpreter is present. It creates a local `.venv`, recreates that `.venv` if it was built with the wrong version or became corrupted, installs the app as a regular package, runs a self-check, and starts MLHeatmap. It does not install Python packages into the user's global Python environment.

After installation, you can start the app again with:

```bash
bash ./run-macos.sh
```

Extra CLI flags still work:

```bash
bash ./run-macos.sh --no-browser
```

Native core install includes:

- Upload and gene mapping
- Normalization
- Interactive heatmaps
- Random Forest
- Logistic Regression (L1)
- Linear SVM
- SHAP analysis
- Excel export
- Browser-side PNG/SVG export

## Docker Install (Full Model Set)

Docker is the recommended path when you want the full model set, including XGBoost and LightGBM.

```bash
docker build -t mlheatmap .
docker run --rm -p 8765:8765 mlheatmap
```

Then open `http://127.0.0.1:8765`.

Docker/full install adds:

- XGBoost
- LightGBM

## Advanced Native Option

If your local environment already supports the optional ML wheels, you can install the full dependency set natively:

```bash
pip install "mlheatmap[full]"
```

## Notes

- Gene mapping tables for human and mouse are packaged with the app.
- Image exports are generated in the browser to avoid Kaleido/Chrome runtime issues.
- Large heatmaps can still use server-side rendering for reliable export and viewing.
