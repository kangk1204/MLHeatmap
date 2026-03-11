"""Group assignment API endpoints."""

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

router = APIRouter(tags=["groups"])


class GroupsRequest(BaseModel):
    session_id: str
    groups: dict


class ExcludeRequest(BaseModel):
    session_id: str
    samples: list


@router.post("/groups")
async def set_groups(request: Request, body: GroupsRequest):
    """Set sample group assignments."""
    session = request.app.state.sessions.get(body.session_id)
    if not session:
        return JSONResponse({"error": "Session not found"}, status_code=404)

    session.groups = body.groups
    return {"groups": session.groups, "n_groups": len(session.groups)}


@router.get("/groups")
async def get_groups(request: Request, session_id: str):
    """Get current group assignments."""
    session = request.app.state.sessions.get(session_id)
    if not session:
        return JSONResponse({"error": "Session not found"}, status_code=404)
    return {"groups": session.groups}


@router.post("/groups/exclude")
async def exclude_samples(request: Request, body: ExcludeRequest):
    """Exclude samples from analysis."""
    session = request.app.state.sessions.get(body.session_id)
    if not session:
        return JSONResponse({"error": "Session not found"}, status_code=404)

    session.excluded_samples = list(set(session.excluded_samples + body.samples))

    # Remove excluded samples from groups
    for group in session.groups:
        session.groups[group] = [
            s for s in session.groups[group] if s not in body.samples
        ]

    return {
        "excluded": session.excluded_samples,
        "remaining": [s for s in session.sample_names if s not in session.excluded_samples],
    }


@router.post("/groups/include")
async def include_samples(request: Request, body: ExcludeRequest):
    """Re-include previously excluded samples."""
    session = request.app.state.sessions.get(body.session_id)
    if not session:
        return JSONResponse({"error": "Session not found"}, status_code=404)

    session.excluded_samples = [s for s in session.excluded_samples if s not in body.samples]
    return {
        "excluded": session.excluded_samples,
        "remaining": [s for s in session.sample_names if s not in session.excluded_samples],
    }
