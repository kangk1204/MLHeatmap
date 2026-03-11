#!/usr/bin/env python3
"""Build gene mapping tables from HGNC and MGI databases.

Downloads and compresses gene mapping data for human and mouse.
Output: src/mlheatmap/data/human_genes.tsv.gz and mouse_genes.tsv.gz
"""

import csv
import gzip
import io
import os
import sys
from pathlib import Path
from urllib.request import urlopen, Request

DATA_DIR = Path(__file__).parent.parent / "src" / "mlheatmap" / "data"


def download(url: str) -> str:
    """Download URL content as string."""
    print(f"  Downloading: {url[:80]}...")
    req = Request(url, headers={"User-Agent": "mlheatmap/0.1"})
    with urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def build_human_genes():
    """Build human gene table from HGNC."""
    print("Building human gene table...")

    # HGNC complete set (TSV)
    url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    try:
        content = download(url)
    except Exception as e:
        print(f"  Failed to download HGNC: {e}")
        print("  Creating minimal human gene table from built-in data...")
        _write_builtin_human()
        return

    reader = csv.DictReader(io.StringIO(content), delimiter="\t")
    output_path = DATA_DIR / "human_genes.tsv.gz"

    count = 0
    with gzip.open(output_path, "wt", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["ensembl_id", "entrez_id", "refseq_ids", "symbol"])

        for row in reader:
            symbol = row.get("symbol", "").strip()
            if not symbol:
                continue

            ensembl = row.get("ensembl_gene_id", "").strip()
            entrez = row.get("entrez_id", "").strip()

            # RefSeq IDs from HGNC
            refseq_accession = row.get("refseq_accession", "").strip()
            refseq_ids = refseq_accession.replace('"', '').strip()

            writer.writerow([ensembl, entrez, refseq_ids, symbol])
            count += 1

    print(f"  Written {count} human genes to {output_path}")


def build_mouse_genes():
    """Build mouse gene table from MGI."""
    print("Building mouse gene table...")

    # MGI gene list with coordinates
    url = "https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt"
    try:
        content = download(url)
    except Exception as e:
        print(f"  Failed to download MGI: {e}")
        print("  Creating minimal mouse gene table from built-in data...")
        _write_builtin_mouse()
        return

    output_path = DATA_DIR / "mouse_genes.tsv.gz"

    count = 0
    with gzip.open(output_path, "wt", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["ensembl_id", "entrez_id", "refseq_ids", "symbol"])

        for line in content.strip().split("\n"):
            if line.startswith("#") or line.startswith("1. MGI"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue

            # MGI format (0-indexed columns):
            #  0: MGI accession id
            #  1: marker type
            #  2: marker symbol
            #  3: marker name
            #  4: genome build
            #  5: Entrez gene id
            #  6: NCBI gene chromosome
            #  7-9: NCBI gene start/end/strand
            # 10: Ensembl gene id
            # 11-14: Ensembl gene chromosome/start/end/strand
            try:
                mgi_type = parts[1].strip() if len(parts) > 1 else ""
                symbol = parts[2].strip() if len(parts) > 2 else ""
                entrez = parts[5].strip() if len(parts) > 5 else ""
                ensembl = parts[10].strip() if len(parts) > 10 else ""

                if not symbol or mgi_type not in ("Gene", "protein coding gene", ""):
                    continue

                writer.writerow([ensembl, entrez, "", symbol])
                count += 1
            except (IndexError, ValueError):
                continue

    print(f"  Written {count} mouse genes to {output_path}")


def _write_builtin_human():
    """Write a minimal human gene table with common genes."""
    output_path = DATA_DIR / "human_genes.tsv.gz"
    # Common human genes for fallback
    genes = [
        ("ENSG00000141510", "7157", "NM_000546", "TP53"),
        ("ENSG00000171862", "7422", "NM_004360", "PTEN"),
        ("ENSG00000141736", "2064", "NM_005228", "ERBB2"),
        ("ENSG00000157764", "673", "NM_004333", "BRAF"),
        ("ENSG00000133703", "3845", "NM_004985", "KRAS"),
        ("ENSG00000149311", "4609", "NM_002467", "ATM"),
        ("ENSG00000012048", "672", "NM_007294", "BRCA1"),
        ("ENSG00000139618", "675", "NM_000059", "BRCA2"),
        ("ENSG00000181449", "4233", "NM_002392", "SOX2"),
        ("ENSG00000111640", "2597", "NM_000601", "GAPDH"),
        ("ENSG00000075624", "60", "NM_001101", "ACTB"),
    ]

    with gzip.open(output_path, "wt", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["ensembl_id", "entrez_id", "refseq_ids", "symbol"])
        for row in genes:
            writer.writerow(row)

    print(f"  Written {len(genes)} common genes (minimal fallback)")


def _write_builtin_mouse():
    """Write a minimal mouse gene table with common genes."""
    output_path = DATA_DIR / "mouse_genes.tsv.gz"
    genes = [
        ("ENSMUSG00000059552", "22059", "", "Trp53"),
        ("ENSMUSG00000013663", "19211", "", "Pten"),
        ("ENSMUSG00000032855", "13866", "", "Erbb2"),
        ("ENSMUSG00000030265", "109880", "", "Braf"),
        ("ENSMUSG00000030265", "16653", "", "Kras"),
        ("ENSMUSG00000029580", "14265", "", "Gapdh"),
        ("ENSMUSG00000029580", "11461", "", "Actb"),
    ]

    with gzip.open(output_path, "wt", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["ensembl_id", "entrez_id", "refseq_ids", "symbol"])
        for row in genes:
            writer.writerow(row)

    print(f"  Written {len(genes)} common genes (minimal fallback)")


if __name__ == "__main__":
    os.makedirs(DATA_DIR, exist_ok=True)
    build_human_genes()
    build_mouse_genes()
    print("\nDone! Gene tables built successfully.")
