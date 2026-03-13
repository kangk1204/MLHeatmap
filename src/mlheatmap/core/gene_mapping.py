"""Gene ID detection and mapping to gene symbols."""

from __future__ import annotations

import csv
import gzip
import importlib.resources
import re
from functools import lru_cache
from typing import Iterable


def detect_id_type(gene_ids: list[str]) -> tuple[str, str]:
    """Detect species and ID type from gene IDs."""
    sample = gene_ids[:200]
    n = len(sample)
    if n == 0:
        return ("unknown", "unknown")

    ensembl_human = sum(1 for gene_id in sample if re.match(r"^ENSG\d{8,}", str(gene_id)))
    ensembl_mouse = sum(1 for gene_id in sample if re.match(r"^ENSMUSG\d{8,}", str(gene_id)))
    refseq = sum(1 for gene_id in sample if re.match(r"^[NX][MR]_\d+", str(gene_id)))
    numeric = sum(1 for gene_id in sample if re.match(r"^\d+$", str(gene_id)))

    threshold = n * 0.4

    if ensembl_human > threshold:
        return ("human", "ensembl")
    if ensembl_mouse > threshold:
        return ("mouse", "ensembl")
    if refseq > threshold:
        return ("unknown", "refseq")
    if numeric > threshold:
        return ("unknown", "entrez")

    alpha_ids = [gene_id for gene_id in sample if re.match(r"^[A-Za-z][A-Za-z0-9\-]*$", str(gene_id))]
    if len(alpha_ids) > threshold:
        uppercase = sum(1 for gene_id in alpha_ids if gene_id == gene_id.upper())
        if uppercase > len(alpha_ids) * 0.5:
            return ("human", "symbol")
        return ("mouse", "symbol")

    return ("unknown", "unknown")


def _gene_table_candidates(species: str) -> Iterable[importlib.resources.abc.Traversable]:
    package_root = importlib.resources.files("mlheatmap.data")
    for extension in (".tsv", ".tsv.gz"):
        yield package_root.joinpath(f"{species}_genes{extension}")


def _resolve_gene_table(species: str):
    for ref in _gene_table_candidates(species):
        if ref.is_file():
            return ref
    return None


def has_gene_table(species: str) -> bool:
    """Return True when the packaged mapping table is available."""
    return _resolve_gene_table(species) is not None


@lru_cache(maxsize=2)
def _load_gene_table(species: str) -> dict[str, dict[str, str]]:
    """Load a packaged gene mapping table. Returns {id_type: {id: symbol}}."""
    ref = _resolve_gene_table(species)
    mapping = {"ensembl": {}, "entrez": {}, "refseq": {}, "symbol": {}}
    if ref is None:
        return mapping

    try:
        with importlib.resources.as_file(ref) as path:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    symbol = row.get("symbol", "").strip()
                    if not symbol:
                        continue

                    ensembl_id = row.get("ensembl_id", "").strip()
                    if ensembl_id:
                        ensembl_base = ensembl_id.split(".")[0]
                        mapping["ensembl"][ensembl_base] = symbol
                        mapping["ensembl"][ensembl_id] = symbol

                    entrez_id = row.get("entrez_id", "").strip()
                    if entrez_id:
                        mapping["entrez"][entrez_id] = symbol

                    for refseq_id in row.get("refseq_ids", "").split("|"):
                        refseq_id = refseq_id.strip()
                        if not refseq_id:
                            continue
                        mapping["refseq"][refseq_id] = symbol
                        mapping["refseq"][refseq_id.split(".")[0]] = symbol

                    mapping["symbol"][symbol] = symbol
                    mapping["symbol"][symbol.upper()] = symbol
        return mapping
    except OSError:
        return {"ensembl": {}, "entrez": {}, "refseq": {}, "symbol": {}}


def map_gene_ids(
    gene_ids: list[str],
    species: str,
    id_type: str = "auto",
) -> tuple[dict[str, str], list[str]]:
    """Map gene IDs to symbols."""
    if species == "unknown":
        for candidate in ["human", "mouse"]:
            mapping, unmapped = map_gene_ids(gene_ids, candidate, id_type)
            if len(unmapped) < len(gene_ids) * 0.5:
                return mapping, unmapped
        return {}, gene_ids

    table = _load_gene_table(species)

    if id_type in {"auto", "unknown"}:
        best_type = "symbol"
        best_count = 0
        for candidate_type in ["ensembl", "entrez", "refseq", "symbol"]:
            count = sum(
                1
                for gene_id in gene_ids[:200]
                if str(gene_id).split(".")[0] in table.get(candidate_type, {})
            )
            if count > best_count:
                best_count = count
                best_type = candidate_type
        id_type = best_type

    lookup = table.get(id_type, {})
    result = {}
    unmapped = []

    for gene_id in gene_ids:
        gene_id_str = str(gene_id).strip()
        gene_id_base = gene_id_str.split(".")[0]
        symbol = lookup.get(gene_id_str) or lookup.get(gene_id_base)
        if symbol:
            result[gene_id_str] = symbol
        else:
            unmapped.append(gene_id_str)

    return result, unmapped
