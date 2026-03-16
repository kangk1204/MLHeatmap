"""In-memory session management."""

from __future__ import annotations

import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd


@dataclass
class Session:
    id: str = field(default_factory=lambda: str(uuid.uuid4()))
    raw_counts: Optional[pd.DataFrame] = None
    mapped_counts: Optional[pd.DataFrame] = None
    normalized: Optional[np.ndarray] = None
    deg_abundance: Optional[np.ndarray] = None
    deg_effect_size_basis: str = ""
    norm_method: str = ""
    size_factors: Optional[np.ndarray] = None
    gene_names: list = field(default_factory=list)
    sample_names: list = field(default_factory=list)
    groups: dict = field(default_factory=dict)
    excluded_samples: list = field(default_factory=list)
    excluded_group_assignments: dict[str, str] = field(default_factory=dict)
    species: str = "unknown"
    id_type: str = "unknown"
    biomarker_results: Optional[dict] = None
    deg_results: Optional[dict] = None
    heatmap_data: Optional[dict] = None
    heatmap_color_scale: str = "RdBu_r"
    metadata: dict[str, Any] = field(default_factory=dict)
    analysis_revision: int = 0
    created_at: float = field(default_factory=time.time)
    last_accessed_at: float = field(default_factory=time.time)
    active_operations: int = 0
    purge_pending: bool = False
    state_lock: threading.RLock = field(default_factory=threading.RLock, repr=False, compare=False)
    _operation_cancel_events: dict[str, threading.Event] = field(
        default_factory=dict,
        repr=False,
        compare=False,
    )

    def touch(self) -> None:
        """Refresh the session access timestamp."""
        with self.state_lock:
            self.last_accessed_at = time.time()

    def revision(self) -> int:
        """Return the current analysis revision."""
        with self.state_lock:
            return self.analysis_revision

    def begin_use(self) -> SessionLease:
        """Acquire an active-operation lease for long-running work."""
        with self.state_lock:
            operation_id = str(uuid.uuid4())
            cancel_event = threading.Event()
            self._operation_cancel_events[operation_id] = cancel_event
            self.active_operations += 1
            self.last_accessed_at = time.time()
            return SessionLease(
                session=self,
                session_id=self.id,
                operation_id=operation_id,
                cancel_event=cancel_event,
            )

    def end_use(self, operation_id: str | None = None) -> tuple[int, bool]:
        """Release an active-operation lease."""
        with self.state_lock:
            if operation_id is not None:
                self._operation_cancel_events.pop(operation_id, None)
            self.active_operations = max(0, self.active_operations - 1)
            self.last_accessed_at = time.time()
            return self.active_operations, self.purge_pending

    def _cancel_active_operations_unlocked(self) -> int:
        for cancel_event in self._operation_cancel_events.values():
            cancel_event.set()
        self.last_accessed_at = time.time()
        return len(self._operation_cancel_events)

    def cancel_active_operations(self) -> int:
        """Signal all active operations for cooperative cancellation."""
        with self.state_lock:
            return self._cancel_active_operations_unlocked()

    def _invalidate_analysis_unlocked(self) -> None:
        """Drop downstream analysis artifacts when the caller already holds the lock."""
        self.analysis_revision += 1
        self.last_accessed_at = time.time()
        self.biomarker_results = None
        self.deg_results = None
        self.heatmap_data = None
        self.metadata.pop("biomarker", None)
        self.metadata.pop("deg", None)

    def _invalidate_normalization_unlocked(self) -> None:
        """Drop normalized data and downstream artifacts when the caller already holds the lock."""
        self.analysis_revision += 1
        self.last_accessed_at = time.time()
        self.normalized = None
        self.deg_abundance = None
        self.deg_effect_size_basis = ""
        self.norm_method = ""
        self.size_factors = None
        self.biomarker_results = None
        self.deg_results = None
        self.heatmap_data = None
        self.metadata.pop("normalization", None)
        self.metadata.pop("biomarker", None)
        self.metadata.pop("deg", None)

    def invalidate_analysis(self) -> None:
        """Drop downstream analysis artifacts when cohort or params change."""
        with self.state_lock:
            self._invalidate_analysis_unlocked()

    def invalidate_normalization(self) -> None:
        """Drop normalized data and all downstream artifacts."""
        with self.state_lock:
            self._invalidate_normalization_unlocked()

    def clear_for_purge(self) -> None:
        """Release heavy session-owned data when the session is purged."""
        with self.state_lock:
            self.raw_counts = None
            self.mapped_counts = None
            self.normalized = None
            self.deg_abundance = None
            self.deg_effect_size_basis = ""
            self.norm_method = ""
            self.size_factors = None
            self.gene_names = []
            self.sample_names = []
            self.groups = {}
            self.excluded_samples = []
            self.excluded_group_assignments = {}
            self.species = "unknown"
            self.id_type = "unknown"
            self.biomarker_results = None
            self.deg_results = None
            self.heatmap_data = None
            self.metadata = {}
            self.purge_pending = True


@dataclass(frozen=True)
class SessionLease:
    session: Session
    session_id: str
    operation_id: str
    cancel_event: threading.Event


class SessionCancelledError(RuntimeError):
    """Raised when a user cancels an active analysis for a session."""


class SessionStore:
    def __init__(self, ttl_hours: int = 4):
        self._sessions: dict[str, Session] = {}
        self._ttl = ttl_hours * 3600
        self._last_cleanup: float = 0.0
        self._lock = threading.RLock()

    def _is_expired(self, session: Session, now: float) -> bool:
        with session.state_lock:
            if session.active_operations > 0:
                return False
            return now - session.last_accessed_at > self._ttl

    def create(self) -> Session:
        self._maybe_cleanup()
        s = Session()
        with self._lock:
            self._sessions[s.id] = s
        return s

    def _finalize_purged_session(self, session_id: str, session: Session) -> None:
        session.clear_for_purge()
        self._sessions.pop(session_id, None)

    def get(self, session_id: str) -> Optional[Session]:
        self._maybe_cleanup()
        now = time.time()
        with self._lock:
            session = self._sessions.get(session_id)
            if session is None:
                return None
            with session.state_lock:
                if session.purge_pending:
                    if session.active_operations == 0:
                        self._finalize_purged_session(session_id, session)
                    return None
            if self._is_expired(session, now):
                del self._sessions[session_id]
                return None
        session.touch()
        return session

    def begin_use(self, session_id: str) -> Optional[SessionLease]:
        """Return a session protected from TTL cleanup during active work."""
        self._maybe_cleanup()
        now = time.time()
        with self._lock:
            session = self._sessions.get(session_id)
            if session is None:
                return None
            with session.state_lock:
                if session.purge_pending:
                    if session.active_operations == 0:
                        self._finalize_purged_session(session_id, session)
                    return None
            if self._is_expired(session, now):
                del self._sessions[session_id]
                return None
            return session.begin_use()

    def end_use(self, session_id: str, operation_id: str | None = None) -> None:
        """Release an active-operation lease if the session still exists."""
        with self._lock:
            session = self._sessions.get(session_id)
        if session is not None:
            remaining, purge_pending = session.end_use(operation_id)
            if purge_pending and remaining == 0:
                with self._lock:
                    current = self._sessions.get(session_id)
                    if current is session:
                        self._finalize_purged_session(session_id, session)

    def cancel(self, session_id: str) -> Optional[dict[str, int | bool]]:
        """Signal all active session work for cooperative cancellation."""
        with self._lock:
            session = self._sessions.get(session_id)
        if session is None:
            return None
        with session.state_lock:
            cancelled = session._cancel_active_operations_unlocked()
            return {
                "cancelled_operations": cancelled,
                "active_operations": session.active_operations,
                "purge_pending": session.purge_pending,
            }

    def purge(self, session_id: str) -> Optional[dict[str, int | bool]]:
        """Mark a session for immediate removal and free it when work stops."""
        with self._lock:
            session = self._sessions.get(session_id)
        if session is None:
            return None
        with session.state_lock:
            session.purge_pending = True
            cancelled = session._cancel_active_operations_unlocked()
            active_operations = session.active_operations
        if active_operations == 0:
            with self._lock:
                current = self._sessions.get(session_id)
                if current is session:
                    self._finalize_purged_session(session_id, session)
        return {
            "purged": active_operations == 0,
            "cancelled_operations": cancelled,
            "active_operations": active_operations,
        }

    def _maybe_cleanup(self):
        now = time.time()
        if now - self._last_cleanup > 300:  # every 5 minutes
            self._cleanup()

    def _cleanup(self):
        now = time.time()
        self._last_cleanup = now
        with self._lock:
            expired = [session_id for session_id, session in self._sessions.items() if self._is_expired(session, now)]
            for session_id in expired:
                del self._sessions[session_id]
