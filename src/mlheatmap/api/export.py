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
    """Export heatmap, plots as images or results as Excel."""
    session = request.app.state.sessions.get(session_id)
    if not session:
        return JSONResponse({"error": "Session not found"}, status_code=404)

    from mlheatmap.core.export import (
        export_auc_image,
        export_heatmap_image,
        export_results_excel,
        export_shap_image,
    )

    try:
        if export_type == "heatmap_png":
            data = export_heatmap_image(session, fmt="png", dpi=dpi)
            return Response(data, media_type="image/png",
                          headers={"Content-Disposition": "attachment; filename=heatmap.png"})
        elif export_type == "heatmap_svg":
            data = export_heatmap_image(session, fmt="svg", dpi=dpi)
            return Response(data, media_type="image/svg+xml",
                          headers={"Content-Disposition": "attachment; filename=heatmap.svg"})
        elif export_type == "shap_png":
            data = export_shap_image(session, fmt="png", dpi=dpi)
            return Response(data, media_type="image/png",
                          headers={"Content-Disposition": "attachment; filename=shap.png"})
        elif export_type == "shap_svg":
            data = export_shap_image(session, fmt="svg", dpi=dpi)
            return Response(data, media_type="image/svg+xml",
                          headers={"Content-Disposition": "attachment; filename=shap.svg"})
        elif export_type == "auc_png":
            data = export_auc_image(session, fmt="png", dpi=dpi)
            return Response(data, media_type="image/png",
                          headers={"Content-Disposition": "attachment; filename=auc.png"})
        elif export_type == "auc_svg":
            data = export_auc_image(session, fmt="svg", dpi=dpi)
            return Response(data, media_type="image/svg+xml",
                          headers={"Content-Disposition": "attachment; filename=auc.svg"})
        elif export_type == "results_excel":
            data = export_results_excel(session)
            return Response(data,
                          media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                          headers={"Content-Disposition": "attachment; filename=mlheatmap_results.xlsx"})
        else:
            return JSONResponse({"error": f"Unknown export type: {export_type}"}, status_code=400)
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=500)
