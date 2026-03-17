"""Shared API validation helpers."""

from __future__ import annotations

import re
import uuid

from fastapi import Request
from fastapi.responses import JSONResponse


INVALID_SESSION_ID_ERROR = "Invalid session_id format. Expected UUID."


def validate_session_id_value(session_id: str) -> tuple[str | None, JSONResponse | None]:
    """Validate and normalize a session UUID string."""
    raw = (session_id or "").strip()
    try:
        normalized = str(uuid.UUID(raw))
    except (AttributeError, TypeError, ValueError):
        return None, JSONResponse({"error": INVALID_SESSION_ID_ERROR}, status_code=400)
    return normalized, None


def get_session_or_error(
    request: Request,
    session_id: str,
    *,
    not_found_message: str = "Session not found",
):
    """Return a validated session or a JSON error response."""
    normalized, error = validate_session_id_value(session_id)
    if error is not None:
        return None, error, None

    session = request.app.state.sessions.get(normalized)
    if session is None:
        return None, JSONResponse({"error": not_found_message}, status_code=404), normalized

    return session, None, normalized


def acquire_session_lease_or_error(
    request: Request,
    session_id: str,
    *,
    not_found_message: str = "Session not found",
):
    """Return a validated session lease with an active-operation cancel token."""
    normalized, error = validate_session_id_value(session_id)
    if error is not None:
        return None, error, None

    lease = request.app.state.sessions.begin_use(normalized)
    if lease is None:
        return None, JSONResponse({"error": not_found_message}, status_code=404), normalized

    return lease, None, normalized


def sanitize_upload_filename(filename: str | None) -> str:
    """Strip path-like components and control characters from upload names."""
    cleaned = (filename or "data.csv").replace("\x00", "").strip()
    cleaned = cleaned.replace("\\", "/").split("/")[-1]
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", cleaned)
    if cleaned in {"", ".", ".."}:
        return "data.csv"
    return cleaned
