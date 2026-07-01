#!/usr/bin/env python3
"""Parse the GSE39582 series matrix + GPL570 platform annotation into a clean,
tumor-only, gene-level expression matrix (samples x genes stays as genes x
samples for R/CMScaller convenience), ready for CMS labeling.

Usage: python scripts/prepare_gse39582.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "generated" / "gse39582" / "sources"
OUT = Path(__file__).resolve().parents[1] / "generated" / "gse39582"


def parse_series_matrix(path: Path):
    source_name = None
    geo_acc = None
    with path.open(encoding="utf-8") as fh:
        header_line = None
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                geo_acc = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                source_name = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                header_line = next(fh)
                break
        assert header_line is not None
        cols = [c.strip('"') for c in header_line.rstrip("\n").split("\t")]
        assert cols[0] == "ID_REF"
        expr = pd.read_csv(fh, sep="\t", header=None, names=cols, index_col=0, na_values=["null"])
        # drop the trailing !series_matrix_table_end line if pandas picked it up as a row
        expr = expr[~expr.index.astype(str).str.startswith("!")]
        expr = expr.apply(pd.to_numeric, errors="coerce")
    meta = pd.DataFrame({"geo_accession": geo_acc, "source_name": source_name}, index=cols[1:])
    return expr, meta


def main() -> int:
    print("Parsing series matrix (large file, may take a minute)...")
    expr, meta = parse_series_matrix(SRC / "GSE39582_series_matrix.txt")
    print(f"  {expr.shape[0]} probes x {expr.shape[1]} samples")

    tumor_mask = meta["source_name"].str.contains("Adenocarcinoma", na=False)
    tumor_samples = meta.index[tumor_mask].tolist()
    print(f"  Tumor samples: {len(tumor_samples)} (of {meta.shape[0]} total)")
    expr = expr[tumor_samples]

    print("Loading GPL570 platform annotation...")
    with (SRC / "GPL570.annot").open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                header = next(fh)
                break
        cols = header.rstrip("\n").split("\t")
        annot = pd.read_csv(fh, sep="\t", header=None, names=cols, dtype=str, on_bad_lines="skip")
        annot = annot[~annot["ID"].astype(str).str.startswith("!")]
    probe_to_symbol = dict(zip(annot["ID"], annot["Gene symbol"]))

    print("Mapping probes to gene symbols (mean across probes per symbol)...")
    symbols = expr.index.to_series().map(probe_to_symbol)
    valid = symbols.notna() & (symbols.astype(str).str.strip() != "") & (~symbols.astype(str).str.contains(r"///"))
    mapped = expr.loc[valid.to_numpy()].copy()
    mapped.index = symbols[valid].to_numpy()
    gene_expr = mapped.groupby(level=0).mean()
    print(f"  {gene_expr.shape[0]} genes x {gene_expr.shape[1]} tumor samples")

    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "gse39582_tumor_gene_expr.tsv.gz"
    gene_expr.to_csv(out_path, sep="\t", compression="gzip")
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
