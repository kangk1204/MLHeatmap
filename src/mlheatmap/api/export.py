"""Export API endpoints."""

from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, Response

router = APIRouter(tags=["export"])


@router.get("/export")
async def export_data(
    request: Request,
    session_id: str = Query(...),
    export_type: str = Query(..., alias="type"),
    dpi: int = Query(300),
):
    """Export workbook data; image exports are handled in the browser."""
    session = request.app.state.sessions.get(session_id)
    if not session:
        return JSONResponse({"error": "Session not found"}, status_code=404)

    from mlheatmap.core.export import export_results_excel

    try:
        if export_type == "results_excel":
            data = export_results_excel(session)
            return Response(
                data,
                media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                headers={"Content-Disposition": "attachment; filename=mlheatmap_results.xlsx"},
            )
        else:
            return JSONResponse(
                {
                    "error": (
                        f"Unsupported export type: {export_type}. "
                        "Image exports are generated in the browser; use results_excel for workbook export."
                    )
                },
                status_code=400,
            )
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=500)
