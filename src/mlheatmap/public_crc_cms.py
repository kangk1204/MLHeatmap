"""Build the public TCGA CRC CMS example used in the manuscript."""

from __future__ import annotations

import argparse
import hashlib
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

# Pinned sha256 of every public source file. The Xena GDC hub and the CRCSC
# label file are mutable URLs, so the cohort is reproducible only against these
# exact snapshots; download_file verifies each file (cached or freshly fetched)
# and fails loudly on any mismatch instead of silently building a different cohort.
EXPECTED_SOURCE_SHA256 = {
    "gdc_coad_star_counts.tsv.gz": "f434486ee869efae683c77ebb6376bd7a11e63ec50e8a62479f73c398db35d42",
    "gdc_read_star_counts.tsv.gz": "48494edbcad75f9478867e069515b3ea0b207dd1f3709cda82e2d6fd4d38f0af",
    "gencode_v36_probemap.tsv": "f027752d663b1da74f4113998bf008e6ebecfa3732607513e027bb6eb1f635b5",
    "crcsc_cms_labels_gold_standard.txt": "9ab613c855f0e31dfdcf417fd5c6f4dfb064cb482132855025e7e3eae7a881b1",
}


def sha256_file(path: Path) -> str:
    """Return the sha256 hex digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksum(spec: DownloadSpec, path: Path) -> None:
    """Raise if a source file's sha256 does not match the pinned snapshot."""
    expected = EXPECTED_SOURCE_SHA256.get(spec.filename)
    if expected is None:
        return
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(
            f"Checksum mismatch for {spec.filename}: expected {expected}, got {actual}. "
            "The upstream public source has changed; the manuscript cohort is pinned to the "
            "recorded sha256 snapshots (see EXPECTED_SOURCE_SHA256)."
        )


def download_file(spec: DownloadSpec, destination: Path, *, force: bool = False) -> Path:
    """Download one public source file unless it is already cached, verifying its sha256."""
    if destination.exists() and not force:
        verify_checksum(spec, destination)
        print(f"  [CACHED] {spec.description}")
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(spec.url, headers={"User-Agent": "MLHeatmap/1.0"})
    print(f"  Downloading {spec.description}...")
    with urllib.request.urlopen(request, timeout=120) as response, temp_path.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    verify_checksum(spec, temp_path)
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
    """Keep one primary-tumor aliquot per TCGA sample barcode, deterministically.

    When a patient has more than one primary-tumor aliquot (e.g. ``01A`` and
    ``01B``), the aliquot is chosen by sorting the full barcodes and taking the
    first (so ``01A`` < ``01B`` < ``01C``), and the output columns are sorted by
    the 15-char sample barcode. Both choices are independent of the source column
    order, so a fresh Xena download that reorders columns yields the same cohort.
    """
    candidates: dict[str, list[str]] = {}
    for column in frame.columns:
        parts = str(column).split("-")
        if len(parts) < 4 or parts[3][:2] != "01":
            continue
        short_barcode = str(column)[:15]
        candidates.setdefault(short_barcode, []).append(str(column))

    short_barcodes = sorted(candidates)
    chosen = [sorted(candidates[key])[0] for key in short_barcodes]
    trimmed = frame.loc[:, chosen].copy()
    trimmed.columns = short_barcodes
    return trimmed


def xena_log2_to_counts(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert Xena log2(count + 1) values back to integer counts."""
    values = np.exp2(frame.to_numpy(dtype=np.float64, copy=False)) - 1.0
    values = np.rint(np.clip(values, a_min=0.0, a_max=None)).astype(np.int64)
    return pd.DataFrame(values, index=frame.index, columns=frame.columns)


def map_unique_gene_symbols(frame: pd.DataFrame, ensembl_to_symbol: dict[str, str]) -> pd.DataFrame:
    """Map Ensembl IDs to gene symbols, summing counts of IDs that share a symbol.

    Several Ensembl IDs can map to one symbol (paralogous/PAR-region entries).
    Their integer counts are summed rather than keeping an arbitrary first row,
    so the per-symbol value is independent of source row order. The output is
    sorted by symbol for a canonical, reproducible gene order.
    """
    symbols = frame.index.to_series().astype(str).map(ensembl_to_symbol.get)
    mask = symbols.notna() & (symbols.astype(str).str.strip() != "")
    mapped = frame.loc[mask.to_numpy()].copy()
    mapped.index = symbols[mask].to_numpy()
    mapped = mapped.groupby(level=0, sort=True).sum()
    return mapped


def apply_low_expression_filter(frame: pd.DataFrame, *, min_count: int = 10, min_fraction: float = 0.2) -> pd.DataFrame:
    """Apply the same low-expression filter used by the manuscript workflow."""
    min_samples = max(2, int(frame.shape[1] * min_fraction))
    mask = (frame >= min_count).sum(axis=1) >= min_samples
    return frame.loc[mask.to_numpy()]


def apply_gold_labels(frame: pd.DataFrame, gold_labels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep labeled samples and rename columns for MLHeatmap group auto-detection."""
    label_counts = gold_labels.groupby("short_barcode")["CMS_gold"].nunique()
    conflicting = sorted(label_counts[label_counts > 1].index.tolist())
    if conflicting:
        raise ValueError(
            f"Conflicting gold CMS labels for {len(conflicting)} barcode(s): {conflicting}. "
            "The gold-label source lists a patient with more than one CMS call; "
            "resolve before building so the assignment is not source-order dependent."
        )
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
    # The low-expression filter is applied across all primary-tumor samples
    # (its >=20%-of-samples denominator is the full primary-tumor column count,
    # not the labeled cohort); gold-label subsetting follows. This order is part
    # of the manuscript workflow and is held fixed for reproducibility.
    filtered = apply_low_expression_filter(mapped)

    print("[6/6] Applying CRCSC gold labels and writing outputs...")
    gold_counts, metadata = apply_gold_labels(filtered, gold_labels)

    # Canonical, source-order-independent layout: samples sorted by id, genes by
    # symbol, so the on-disk order (and thus the downstream CV fold composition)
    # does not depend on how Xena happens to order its columns/rows.
    gold_counts = gold_counts.sort_index(axis=0).sort_index(axis=1)
    metadata = metadata.sort_values("sample_id").reset_index(drop=True)

    cms_distribution = metadata["CMS_gold"].value_counts().sort_index().to_dict()
    expected_distribution = {"CMS1": 76, "CMS2": 220, "CMS3": 72, "CMS4": 143}
    n_samples = int(gold_counts.shape[1])
    if n_samples != 511 or cms_distribution != expected_distribution:
        raise ValueError(
            f"Unexpected cohort: {n_samples} samples, distribution {cms_distribution} "
            f"(expected 511 samples, {expected_distribution}). Source coverage may have drifted."
        )

    counts_path = output_dir / "tcga_crc_cms_gold_counts.tsv.gz"
    metadata_path = output_dir / "tcga_crc_cms_gold_metadata.tsv"
    summary_path = output_dir / "tcga_crc_cms_gold_labels.json"

    # mtime=0 makes the gzip byte-deterministic so the cohort sha256 is a stable anchor.
    gold_counts.to_csv(counts_path, sep="\t", compression={"method": "gzip", "mtime": 0})
    metadata.to_csv(metadata_path, sep="\t", index=False)

    summary = {
        "dataset": "TCGA COAD + READ",
        "source": "UCSC Xena GDC hub",
        "xena_value_encoding": "log2(count + 1)",
        "reverse_transform": "round(2^x - 1)",
        "gold_label_source": GOLD_LABELS_URL,
        "gold_label_column": GOLD_LABELS_COLUMN,
        "n_samples": n_samples,
        "n_genes": int(gold_counts.shape[0]),
        "cms_distribution": {key: int(value) for key, value in cms_distribution.items()},
        "source_sha256": dict(EXPECTED_SOURCE_SHA256),
        "cohort_sha256": sha256_file(counts_path),
        "output_counts": str(counts_path),
        "output_metadata": str(metadata_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"  Counts: {counts_path}")
    print(f"  Metadata: {metadata_path}")
    print(f"  Summary: {summary_path}")
    print(f"  Cohort sha256: {summary['cohort_sha256']}")
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
