"""Normalization API endpoints."""

import asyncio
import numpy as np
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
    from mlheatmap.api.validation import acquire_session_lease_or_error
    from mlheatmap.core.cancellation import raise_if_cancelled
    from mlheatmap.core.normalization import (
        deseq2_normalize,
        log2_normalize,
        tpm_abundance,
        tpm_normalize,
    )

    lease, error, _ = acquire_session_lease_or_error(
        request,
        body.session_id,
        not_found_message="Session not found or no data uploaded",
    )
    if error is not None:
        return error
    session = lease.session

    try:
        with session.state_lock:
            if session.raw_counts is None:
                return JSONResponse({"error": "Session not found or no data uploaded"}, status_code=404)
            source_df = session.mapped_counts if session.mapped_counts is not None else session.raw_counts
            counts = source_df.to_numpy(dtype=np.float64, copy=True)
            gene_names = source_df.index.tolist()
            sample_names = source_df.columns.tolist()

        cancel_check = lease.cancel_event.is_set

        def _run_normalization():
            raise_if_cancelled(cancel_check)
            if body.method == "deseq2":
                normalized, size_factors = deseq2_normalize(counts, return_size_factors=True)
                deg_abundance = counts / size_factors[np.newaxis, :]
                deg_effect_size_basis = "size_factor_normalized_counts"
            elif body.method == "tpm":
                normalized = tpm_normalize(counts)
                size_factors = None
                deg_abundance = tpm_abundance(counts)
                deg_effect_size_basis = "tpm"
            elif body.method == "log2":
                normalized = log2_normalize(counts)
                size_factors = None
                deg_abundance = counts.copy()
                deg_effect_size_basis = "counts"
            else:
                raise ValueError(f"Unknown method: {body.method}")
            raise_if_cancelled(cancel_check)
            return normalized, size_factors, deg_abundance, deg_effect_size_basis

        try:
            normalized, size_factors, deg_abundance, deg_effect_size_basis = await asyncio.to_thread(_run_normalization)
        except ValueError as exc:
            return JSONResponse({"error": str(exc)}, status_code=400)
        except Exception as exc:
            from mlheatmap.api.session import SessionCancelledError

            if isinstance(exc, SessionCancelledError):
                return JSONResponse({"error": str(exc)}, status_code=409)
            raise

        with session.state_lock:
            session._invalidate_analysis_unlocked()
            session.normalized = normalized
            session.size_factors = size_factors
            session.deg_abundance = deg_abundance
            session.deg_effect_size_basis = deg_effect_size_basis
            session.norm_method = body.method
            session.gene_names = gene_names
            session.sample_names = sample_names
            session.metadata["normalization"] = {
                "method": body.method,
                "shape": [int(normalized.shape[0]), int(normalized.shape[1])],
                "effect_size_basis": session.deg_effect_size_basis,
            }

        finite = normalized[np.isfinite(normalized)]
        stats = {
            "min": float(np.min(finite)) if len(finite) > 0 else 0,
            "max": float(np.max(finite)) if len(finite) > 0 else 0,
            "median": float(np.median(finite)) if len(finite) > 0 else 0,
            "q25": float(np.percentile(finite, 25)) if len(finite) > 0 else 0,
            "q75": float(np.percentile(finite, 75)) if len(finite) > 0 else 0,
        }

        sample_size = min(1000, len(finite))
        rng = np.random.default_rng(42)
        dist_sample = rng.choice(finite, size=sample_size, replace=False).tolist() if sample_size > 0 else []

        return {
            "method": body.method,
            "shape": list(normalized.shape),
            "stats": stats,
            "distribution_sample": dist_sample,
        }
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
