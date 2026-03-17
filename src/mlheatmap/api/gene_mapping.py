"""Gene mapping API endpoint."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from mlheatmap.api.validation import acquire_session_lease_or_error

router = APIRouter(tags=["gene-mapping"])


class MappingRequest(BaseModel):
    session_id: str
    species: str = "human"
    id_type: str = "auto"


@router.post("/gene-mapping")
async def map_genes(request: Request, body: MappingRequest):
    """Map gene IDs to gene symbols."""
    from mlheatmap.core.gene_mapping import map_gene_ids

    lease, error, _ = acquire_session_lease_or_error(
        request,
        body.session_id,
        not_found_message="Session not found or no data",
    )
    if error is not None:
        return error
    session = lease.session

    try:
        with session.state_lock:
            if session.raw_counts is None:
                return JSONResponse({"error": "Session not found or no data"}, status_code=404)
            raw_counts = session.raw_counts.copy()
            detected_id_type = session.id_type

        gene_ids = raw_counts.index.astype(str).tolist()
        mapping, unmapped = map_gene_ids(gene_ids, body.species, body.id_type)
        warning = None
        preserve_raw_symbols = body.id_type == "symbol" or (body.id_type == "auto" and detected_id_type == "symbol")

        if mapping:
            df = raw_counts.copy()
            df.index = df.index.astype(str)
            if preserve_raw_symbols:
                if unmapped:
                    warning = (
                        "Some uploaded gene symbols were not found in the packaged mapping table. "
                        "They were retained as uploaded identifiers."
                    )
                df["_symbol"] = df.index.map(lambda gene_id: mapping.get(str(gene_id), str(gene_id)))
            else:
                df["_symbol"] = df.index.map(lambda gene_id: mapping.get(str(gene_id)))
                df = df.dropna(subset=["_symbol"])
            df = df.groupby("_symbol").sum()
            df.index.name = None
            mapped_counts = df
        else:
            mapped_counts = raw_counts.copy()
            warning = "No gene mappings were found. Continuing with raw uploaded identifiers."

        with session.state_lock:
            session.mapped_counts = mapped_counts
            session.gene_names = mapped_counts.index.astype(str).tolist()
            session.sample_names = mapped_counts.columns.tolist()
            session.species = body.species
            session._invalidate_normalization_unlocked()
            session.metadata["mapping"] = {
                "species": body.species,
                "id_type": body.id_type,
                "mapped_count": len(mapping),
                "unmapped_count": len(unmapped),
                "total": len(gene_ids),
                "used_raw_ids": not bool(mapping),
            }

        response = {
            "mapped_count": len(mapping),
            "unmapped_count": len(unmapped),
            "total": len(gene_ids),
            "unmapped_sample": unmapped[:50],
        }
        if warning:
            response["warning"] = warning
        return response
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
