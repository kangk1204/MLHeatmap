"""Runtime capability endpoint."""

from fastapi import APIRouter

from mlheatmap.core.capabilities import get_capabilities

router = APIRouter(tags=["capabilities"])


@router.get("/capabilities")
async def capabilities():
    """Return install-time and runtime capabilities for the current app."""
    return get_capabilities()
