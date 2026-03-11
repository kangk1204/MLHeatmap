"""Gene mapping API endpoint."""

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

import pandas as pd

router = APIRouter(tags=["gene-mapping"])


class MappingRequest(BaseModel):
    session_id: str
    species: str = "human"
    id_type: str = "auto"


@router.post("/gene-mapping")
async def map_genes(request: Request, body: MappingRequest):
    """Map gene IDs to gene symbols."""
    from mlheatmap.core.gene_mapping import map_gene_ids

    session = request.app.state.sessions.get(body.session_id)
    if not session or session.raw_counts is None:
        return JSONResponse({"error": "Session not found or no data"}, status_code=404)

    gene_ids = session.raw_counts.index.astype(str).tolist()
    mapping, unmapped = map_gene_ids(gene_ids, body.species, body.id_type)

    if mapping:
        # Remap the DataFrame
        df = session.raw_counts.copy()
        df.index = df.index.astype(str)

        # Map gene IDs to symbols
        df["_symbol"] = df.index.map(lambda x: mapping.get(str(x), None))
        df = df.dropna(subset=["_symbol"])

        # Aggregate duplicates by summing
        df = df.groupby("_symbol").sum()
        df.index.name = None

        session.mapped_counts = df
        session.gene_names = df.index.tolist()
        session.sample_names = df.columns.tolist()
        session.species = body.species
    else:
        # No mapping available, use raw IDs
        session.mapped_counts = session.raw_counts.copy()
        session.gene_names = session.raw_counts.index.astype(str).tolist()

    return {
        "mapped_count": len(mapping),
        "unmapped_count": len(unmapped),
        "total": len(gene_ids),
        "unmapped_sample": unmapped[:50],
    }
