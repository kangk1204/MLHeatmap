"""Gene ID detection and mapping to gene symbols."""

import gzip
import importlib.resources
import re
from functools import lru_cache
from typing import Optional


def detect_id_type(gene_ids: list[str]) -> tuple[str, str]:
    """Detect species and ID type from gene IDs.

    Returns (species, id_type) where:
      species: "human", "mouse", or "unknown"
      id_type: "ensembl", "entrez", "refseq", "symbol", or "unknown"
    """
    sample = gene_ids[:200]
    n = len(sample)
    if n == 0:
        return ("unknown", "unknown")

    ensembl_human = sum(1 for g in sample if re.match(r"^ENSG\d{8,}", str(g)))
    ensembl_mouse = sum(1 for g in sample if re.match(r"^ENSMUSG\d{8,}", str(g)))
    refseq = sum(1 for g in sample if re.match(r"^[NX][MR]_\d+", str(g)))
    numeric = sum(1 for g in sample if re.match(r"^\d+$", str(g)))

    threshold = n * 0.4

    if ensembl_human > threshold:
        return ("human", "ensembl")
    if ensembl_mouse > threshold:
        return ("mouse", "ensembl")
    if refseq > threshold:
        return ("unknown", "refseq")
    if numeric > threshold:
        return ("unknown", "entrez")

    # Check if already gene symbols
    alpha_ids = [g for g in sample if re.match(r"^[A-Za-z][A-Za-z0-9\-]*$", str(g))]
    if len(alpha_ids) > threshold:
        uppercase = sum(1 for g in alpha_ids if g == g.upper())
        if uppercase > len(alpha_ids) * 0.5:
            return ("human", "symbol")
        return ("mouse", "symbol")

    return ("unknown", "unknown")


@lru_cache(maxsize=2)
def _load_gene_table(species: str) -> dict[str, dict[str, str]]:
    """Load compressed gene mapping table. Returns {id_type: {id: symbol}}."""
    filename = f"{species}_genes.tsv.gz"

    try:
        ref = importlib.resources.files("mlheatmap.data").joinpath(filename)
        mapping = {"ensembl": {}, "entrez": {}, "refseq": {}, "symbol": {}}

        with importlib.resources.as_file(ref) as path:
            with gzip.open(path, "rt") as f:
                header = f.readline().strip().split("\t")
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 4:
                        continue
                    row = dict(zip(header, parts))
                    symbol = row.get("symbol", "").strip()
                    if not symbol:
                        continue

                    ens = row.get("ensembl_id", "").strip()
                    if ens:
                        # Strip version
                        ens_base = ens.split(".")[0]
                        mapping["ensembl"][ens_base] = symbol
                        mapping["ensembl"][ens] = symbol

                    entrez = row.get("entrez_id", "").strip()
                    if entrez:
                        mapping["entrez"][entrez] = symbol

                    for rid in row.get("refseq_ids", "").split("|"):
                        rid = rid.strip()
                        if rid:
                            mapping["refseq"][rid] = symbol
                            mapping["refseq"][rid.split(".")[0]] = symbol

                    # Also map symbol to itself (for already-symbol inputs)
                    mapping["symbol"][symbol] = symbol
                    mapping["symbol"][symbol.upper()] = symbol

        return mapping
    except Exception:
        return {"ensembl": {}, "entrez": {}, "refseq": {}, "symbol": {}}


def map_gene_ids(
    gene_ids: list[str],
    species: str,
    id_type: str = "auto",
) -> tuple[dict[str, str], list[str]]:
    """Map gene IDs to symbols.

    Returns (mapping_dict {original_id: symbol}, list_of_unmapped_ids).
    """
    if species == "unknown":
        # Try both
        for sp in ["human", "mouse"]:
            mapping, unmapped = map_gene_ids(gene_ids, sp, id_type)
            if len(unmapped) < len(gene_ids) * 0.5:
                return mapping, unmapped
        return {}, gene_ids

    table = _load_gene_table(species)

    if id_type == "auto" or id_type == "unknown":
        # Try each type, pick best
        best_type = "symbol"
        best_count = 0
        for t in ["ensembl", "entrez", "refseq", "symbol"]:
            count = sum(1 for g in gene_ids[:200] if str(g).split(".")[0] in table.get(t, {}))
            if count > best_count:
                best_count = count
                best_type = t
        id_type = best_type

    lookup = table.get(id_type, {})
    result = {}
    unmapped = []

    for gid in gene_ids:
        gid_str = str(gid).strip()
        gid_base = gid_str.split(".")[0]

        symbol = lookup.get(gid_str) or lookup.get(gid_base)
        if symbol:
            result[gid_str] = symbol
        else:
            unmapped.append(gid_str)

    return result, unmapped
