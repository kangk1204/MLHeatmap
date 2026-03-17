"""Build the public TCGA CRC CMS example used in the manuscript."""

from __future__ import annotations

import argparse
import json
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


COAD_COUNTS_URL = "https://gdc.xenahubs.net/download/TCGA-COAD.star_counts.tsv.gz"
READ_COUNTS_URL = "https://gdc.xenahubs.net/download/TCGA-READ.star_counts.tsv.gz"
PROBEMAP_URL = "https://gdc.xenahubs.net/download/gencode.v36.annotation.gtf.gene.probemap"
GOLD_LABELS_URL = (
    "https://raw.githubusercontent.com/Sage-Bionetworks/crc-cms-kras/"
    "master/020717/cms_labels_public_all.txt"
)
GOLD_LABELS_COLUMN = "CMS_final_network_plus_RFclassifier_in_nonconsensus_samples"
SUMMARY_PREFIXES = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous", "__")


@dataclass(frozen=True)
class DownloadSpec:
    key: str
    filename: str
    url: str
    description: str


DOWNLOAD_SPECS = (
    DownloadSpec("coad_counts", "gdc_coad_star_counts.tsv.gz", COAD_COUNTS_URL, "TCGA COAD STAR counts"),
    DownloadSpec("read_counts", "gdc_read_star_counts.tsv.gz", READ_COUNTS_URL, "TCGA READ STAR counts"),
    DownloadSpec("probemap", "gencode_v36_probemap.tsv", PROBEMAP_URL, "Gencode v36 probe map"),
    DownloadSpec("gold_labels", "crcsc_cms_labels_gold_standard.txt", GOLD_LABELS_URL, "CRCSC gold CMS labels"),
)


def download_file(spec: DownloadSpec, destination: Path, *, force: bool = False) -> Path:
    """Download one public source file unless it is already cached."""
    if destination.exists() and not force:
        print(f"  [CACHED] {spec.description}")
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(spec.url, headers={"User-Agent": "MLHeatmap/1.0"})
    print(f"  Downloading {spec.description}...")
    with urllib.request.urlopen(request, timeout=120) as response, temp_path.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    temp_path.replace(destination)
    size_mb = destination.stat().st_size / (1024 * 1024)
    print(f"  [OK] {destination.name} ({size_mb:.1f} MB)")
    return destination


def download_sources(download_dir: Path, *, force: bool = False) -> dict[str, Path]:
    """Download all public source files and return their local paths."""
    paths: dict[str, Path] = {}
    for spec in DOWNLOAD_SPECS:
        paths[spec.key] = download_file(spec, download_dir / spec.filename, force=force)
    return paths


def load_probemap(path: Path) -> dict[str, str]:
    """Load the Xena probe map as an Ensembl-to-symbol lookup."""
    frame = pd.read_csv(path, sep="\t", usecols=["id", "gene"], dtype=str)
    frame = frame.dropna(subset=["id", "gene"])
    frame["id"] = frame["id"].str.strip()
    frame["gene"] = frame["gene"].str.strip()
    frame = frame[(frame["id"] != "") & (frame["gene"] != "")]
    return dict(zip(frame["id"], frame["gene"]))


def load_gold_labels(path: Path) -> pd.DataFrame:
    """Load CRCSC final CMS labels for TCGA and exclude NOLBL samples."""
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if GOLD_LABELS_COLUMN not in frame.columns:
        raise ValueError(f"Missing gold-label column: {GOLD_LABELS_COLUMN}")

    dataset = frame["dataset"].str.lower() if "dataset" in frame.columns else pd.Series("", index=frame.index)
    tcga = frame[dataset == "tcga"].copy()
    tcga = tcga.rename(columns={GOLD_LABELS_COLUMN: "CMS_gold"})
    tcga = tcga[tcga["CMS_gold"].notna()].copy()
    tcga = tcga[tcga["CMS_gold"] != "NOLBL"].copy()
    tcga["short_barcode"] = tcga["sample"].str.slice(0, 12)
    return tcga[["short_barcode", "CMS_gold"]]


def load_xena_counts(path: Path) -> pd.DataFrame:
    """Load one Xena STAR matrix."""
    return pd.read_csv(path, sep="\t", index_col=0, low_memory=False)


def remove_summary_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop STAR summary rows and keep gene features only."""
    mask = ~frame.index.to_series().astype(str).str.startswith(SUMMARY_PREFIXES)
    return frame.loc[mask.to_numpy()]


def select_primary_tumor_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep one primary-tumor aliquot per TCGA sample barcode."""
    selected: dict[str, str] = {}
    for column in frame.columns:
        parts = str(column).split("-")
        if len(parts) < 4 or parts[3][:2] != "01":
            continue
        short_barcode = str(column)[:15]
        selected.setdefault(short_barcode, str(column))

    trimmed = frame.loc[:, list(selected.values())].copy()
    trimmed.columns = list(selected.keys())
    return trimmed


def xena_log2_to_counts(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert Xena log2(count + 1) values back to integer counts."""
    values = np.exp2(frame.to_numpy(dtype=np.float64, copy=False)) - 1.0
    values = np.rint(np.clip(values, a_min=0.0, a_max=None)).astype(np.int64)
    return pd.DataFrame(values, index=frame.index, columns=frame.columns)


def map_unique_gene_symbols(frame: pd.DataFrame, ensembl_to_symbol: dict[str, str]) -> pd.DataFrame:
    """Map Ensembl IDs to gene symbols and keep the first row per unique symbol."""
    kept_rows: list[str] = []
    symbols: list[str] = []
    seen: set[str] = set()
    for ensembl_id in frame.index:
        symbol = ensembl_to_symbol.get(str(ensembl_id))
        if not symbol or symbol in seen:
            continue
        kept_rows.append(str(ensembl_id))
        symbols.append(symbol)
        seen.add(symbol)

    mapped = frame.loc[kept_rows].copy()
    mapped.index = symbols
    return mapped


def apply_low_expression_filter(frame: pd.DataFrame, *, min_count: int = 10, min_fraction: float = 0.2) -> pd.DataFrame:
    """Apply the same low-expression filter used by the manuscript workflow."""
    min_samples = max(2, int(frame.shape[1] * min_fraction))
    mask = (frame >= min_count).sum(axis=1) >= min_samples
    return frame.loc[mask.to_numpy()]


def apply_gold_labels(frame: pd.DataFrame, gold_labels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep labeled samples and rename columns for MLHeatmap group auto-detection."""
    label_map = dict(zip(gold_labels["short_barcode"], gold_labels["CMS_gold"]))
    selected_columns: list[str] = []
    renamed_columns: dict[str, str] = {}
    metadata_rows: list[dict[str, str]] = []

    for barcode in frame.columns:
        short_barcode = str(barcode)[:12]
        cms_label = label_map.get(short_barcode)
        if cms_label is None:
            continue
        sample_id = f"{cms_label}__{str(barcode).replace('TCGA-', '').replace('-', '')[:8]}"
        selected_columns.append(str(barcode))
        renamed_columns[str(barcode)] = sample_id
        metadata_rows.append(
            {
                "tcga_barcode": str(barcode),
                "short_barcode": short_barcode,
                "CMS_gold": cms_label,
                "sample_id": sample_id,
            }
        )

    labeled = frame.loc[:, selected_columns].rename(columns=renamed_columns)
    metadata = pd.DataFrame(metadata_rows, columns=["tcga_barcode", "short_barcode", "CMS_gold", "sample_id"])
    return labeled, metadata


def build_public_crc_cms_example(output_dir: Path, *, force_download: bool = False) -> dict[str, object]:
    """Download public sources and create the manuscript-ready TCGA CRC CMS matrix."""
    output_dir.mkdir(parents=True, exist_ok=True)
    download_dir = output_dir / "sources"

    print("[1/6] Downloading public sources...")
    paths = download_sources(download_dir, force=force_download)

    print("[2/6] Loading public metadata...")
    ensembl_to_symbol = load_probemap(paths["probemap"])
    gold_labels = load_gold_labels(paths["gold_labels"])

    print("[3/6] Loading TCGA Xena STAR matrices...")
    coad = load_xena_counts(paths["coad_counts"])
    read = load_xena_counts(paths["read_counts"])
    merged = pd.concat([coad, read], axis=1, join="inner")

    print("[4/6] Filtering TCGA primary tumors and reversing Xena transform...")
    merged = remove_summary_rows(merged)
    merged = select_primary_tumor_columns(merged)
    raw_counts = xena_log2_to_counts(merged)

    print("[5/6] Mapping genes and applying manuscript filter...")
    mapped = map_unique_gene_symbols(raw_counts, ensembl_to_symbol)
    filtered = apply_low_expression_filter(mapped)

    print("[6/6] Applying CRCSC gold labels and writing outputs...")
    gold_counts, metadata = apply_gold_labels(filtered, gold_labels)
    counts_path = output_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    metadata_path = output_dir / "tcga_crc_cms_gold_metadata.tsv"
    summary_path = output_dir / "tcga_crc_cms_gold_labels.json"

    gold_counts.to_csv(counts_path, sep="\t", compression="gzip")
    metadata.to_csv(metadata_path, sep="\t", index=False)

    cms_distribution = metadata["CMS_gold"].value_counts().sort_index().to_dict()
    summary = {
        "dataset": "TCGA COAD + READ",
        "source": "UCSC Xena GDC hub",
        "xena_value_encoding": "log2(count + 1)",
        "reverse_transform": "round(2^x - 1)",
        "gold_label_source": GOLD_LABELS_URL,
        "gold_label_column": GOLD_LABELS_COLUMN,
        "n_samples": int(gold_counts.shape[1]),
        "n_genes": int(gold_counts.shape[0]),
        "cms_distribution": {key: int(value) for key, value in cms_distribution.items()},
        "output_counts": str(counts_path),
        "output_metadata": str(metadata_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"  Counts: {counts_path}")
    print(f"  Metadata: {metadata_path}")
    print(f"  Summary: {summary_path}")
    print(f"  Final matrix: {gold_counts.shape[0]} genes x {gold_counts.shape[1]} samples")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="mlheatmap-download-crc-cms",
        description="Download and build the public TCGA CRC CMS gold-label example for MLHeatmap.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("generated") / "crc_cms_public",
        help="Directory where downloads and processed outputs will be written.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download public source files even if cached.",
    )
    args = parser.parse_args(argv)

    build_public_crc_cms_example(args.output_dir, force_download=args.force_download)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
