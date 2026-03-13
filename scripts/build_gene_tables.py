#!/usr/bin/env python3
"""Build packaged gene mapping tables from HGNC and MGI databases.

Downloads and writes:
  - src/mlheatmap/data/human_genes.tsv
  - src/mlheatmap/data/mouse_genes.tsv
"""

from __future__ import annotations

import csv
import io
import os
from pathlib import Path
from urllib.request import Request, urlopen

DATA_DIR = Path(__file__).parent.parent / "src" / "mlheatmap" / "data"


def download(url: str) -> str:
    """Download URL content as text."""
    print(f"Downloading: {url[:80]}...")
    request = Request(url, headers={"User-Agent": "mlheatmap/0.1"})
    with urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def build_human_genes() -> None:
    """Build the human gene table from HGNC."""
    print("Building human gene table...")
    content = download(
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    )
    reader = csv.DictReader(io.StringIO(content), delimiter="\t")
    output_path = DATA_DIR / "human_genes.tsv"

    count = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ensembl_id", "entrez_id", "refseq_ids", "symbol"])
        for row in reader:
            symbol = row.get("symbol", "").strip()
            if not symbol:
                continue
            writer.writerow(
                [
                    row.get("ensembl_gene_id", "").strip(),
                    row.get("entrez_id", "").strip(),
                    row.get("refseq_accession", "").replace('"', "").strip(),
                    symbol,
                ]
            )
            count += 1

    print(f"Wrote {count} human genes to {output_path}")


def build_mouse_genes() -> None:
    """Build the mouse gene table from MGI."""
    print("Building mouse gene table...")
    content = download("https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt")
    output_path = DATA_DIR / "mouse_genes.tsv"

    count = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ensembl_id", "entrez_id", "refseq_ids", "symbol"])
        for line in content.splitlines():
            if line.startswith("#") or line.startswith("1. MGI"):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue

            marker_type = parts[1].strip()
            symbol = parts[2].strip()
            entrez = parts[5].strip()
            ensembl = parts[10].strip()
            if not symbol or marker_type not in ("Gene", "protein coding gene", ""):
                continue

            writer.writerow([ensembl, entrez, "", symbol])
            count += 1

    print(f"Wrote {count} mouse genes to {output_path}")


if __name__ == "__main__":
    os.makedirs(DATA_DIR, exist_ok=True)
    build_human_genes()
    build_mouse_genes()
    print("Done.")
