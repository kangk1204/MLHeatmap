"""Session lifecycle control endpoints."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from mlheatmap.api.validation import validate_session_id_value

router = APIRouter(tags=["session"])


class SessionActionRequest(BaseModel):
    session_id: str


@router.post("/session/cancel")
async def cancel_session_work(request: Request, body: SessionActionRequest):
    """Signal cooperative cancellation for active work tied to a session."""
    normalized_session_id, error = validate_session_id_value(body.session_id)
    if error is not None:
        return error

    result = request.app.state.sessions.cancel(normalized_session_id)
    if result is None:
        return JSONResponse({"error": "Session not found"}, status_code=404)
    return result


@router.post("/session/purge")
async def purge_session(request: Request, body: SessionActionRequest):
    """Immediately retire a session and free it once any active work stops."""
    normalized_session_id, error = validate_session_id_value(body.session_id)
    if error is not None:
        return error

    result = request.app.state.sessions.purge(normalized_session_id)
    if result is None:
        return JSONResponse({"error": "Session not found"}, status_code=404)
    return result
