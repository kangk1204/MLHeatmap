"""Export API endpoints."""

import copy
from types import SimpleNamespace

import numpy as np
from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, Response

from mlheatmap.api.validation import get_session_or_error

router = APIRouter(tags=["export"])


def _snapshot_export_session(session):
    """Copy the session state needed for a consistent workbook export."""
    with session.state_lock:
        return SimpleNamespace(
            id=session.id,
            species=session.species,
            id_type=session.id_type,
            norm_method=session.norm_method,
            deg_effect_size_basis=session.deg_effect_size_basis,
            gene_names=list(session.gene_names),
            sample_names=list(session.sample_names),
            excluded_samples=list(session.excluded_samples),
            groups=copy.deepcopy(session.groups),
            metadata=copy.deepcopy(session.metadata),
            normalized=None if session.normalized is None else np.array(session.normalized, copy=True),
            biomarker_results=copy.deepcopy(session.biomarker_results),
            deg_results=copy.deepcopy(session.deg_results),
        )


@router.get("/export")
async def export_data(
    request: Request,
    session_id: str = Query(...),
    export_type: str = Query(..., alias="type"),
    dpi: int = Query(300),
):
    """Export workbook data; image exports are handled in the browser."""
    session, error, _ = get_session_or_error(request, session_id)
    if error is not None:
        return error

    from mlheatmap.core.export import export_results_excel

    try:
        if export_type == "results_excel":
            data = export_results_excel(_snapshot_export_session(session), include_normalized_expression=False)
            return Response(
                data,
                media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                headers={"Content-Disposition": "attachment; filename=mlheatmap_results_compact.xlsx"},
            )
        if export_type == "results_excel_full":
            data = export_results_excel(_snapshot_export_session(session), include_normalized_expression=True)
            return Response(
                data,
                media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                headers={"Content-Disposition": "attachment; filename=mlheatmap_results_full.xlsx"},
            )
        return JSONResponse(
            {
                "error": (
                    f"Unsupported export type: {export_type}. "
                    "Image exports are generated in the browser; use results_excel or results_excel_full for workbook export."
                )
            },
            status_code=400,
        )
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)
