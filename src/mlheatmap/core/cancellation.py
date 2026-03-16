"""Helpers for cooperative cancellation of long-running analyses."""

from __future__ import annotations

from collections.abc import Callable

from mlheatmap.api.session import SessionCancelledError


def raise_if_cancelled(cancel_check: Callable[[], bool] | None) -> None:
    """Raise a shared cancellation error when the caller requests termination."""
    if cancel_check and cancel_check():
        raise SessionCancelledError("Analysis cancelled at user request.")
