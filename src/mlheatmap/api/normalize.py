"""Normalization API endpoints."""

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

router = APIRouter(tags=["normalize"])


class NormalizeRequest(BaseModel):
    session_id: str
    method: str = "log2"


@router.post("/normalize")
async def normalize(request: Request, body: NormalizeRequest):
    """Normalize count matrix."""
    import numpy as np
    from mlheatmap.core.normalization import deseq2_normalize, log2_normalize, tpm_normalize

    session = request.app.state.sessions.get(body.session_id)
    if not session or session.raw_counts is None:
        return JSONResponse({"error": "Session not found or no data uploaded"}, status_code=404)

    df = session.mapped_counts if session.mapped_counts is not None else session.raw_counts
    counts = df.values.astype(float)

    if body.method == "deseq2":
        normalized = deseq2_normalize(counts)
    elif body.method == "tpm":
        normalized = tpm_normalize(counts)
    elif body.method == "log2":
        normalized = log2_normalize(counts)
    else:
        return JSONResponse({"error": f"Unknown method: {body.method}"}, status_code=400)

    session.normalized = normalized
    session.norm_method = body.method
    session.gene_names = df.index.tolist()
    session.sample_names = df.columns.tolist()

    finite = normalized[np.isfinite(normalized)]
    stats = {
        "min": float(np.min(finite)) if len(finite) > 0 else 0,
        "max": float(np.max(finite)) if len(finite) > 0 else 0,
        "median": float(np.median(finite)) if len(finite) > 0 else 0,
        "q25": float(np.percentile(finite, 25)) if len(finite) > 0 else 0,
        "q75": float(np.percentile(finite, 75)) if len(finite) > 0 else 0,
    }

    # Distribution data for plotting (sample 1000 values)
    sample_size = min(1000, len(finite))
    rng = np.random.default_rng(42)
    dist_sample = rng.choice(finite, size=sample_size, replace=False).tolist()

    return {
        "method": body.method,
        "shape": list(normalized.shape),
        "stats": stats,
        "distribution_sample": dist_sample,
    }
