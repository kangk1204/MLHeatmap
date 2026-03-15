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
    state_lock: threading.RLock = field(default_factory=threading.RLock, repr=False, compare=False)

    def touch(self) -> None:
        """Refresh the session access timestamp."""
        with self.state_lock:
            self.last_accessed_at = time.time()

    def revision(self) -> int:
        """Return the current analysis revision."""
        with self.state_lock:
            return self.analysis_revision

    def begin_use(self) -> None:
        """Acquire an active-operation lease for long-running work."""
        with self.state_lock:
            self.active_operations += 1
            self.last_accessed_at = time.time()

    def end_use(self) -> None:
        """Release an active-operation lease."""
        with self.state_lock:
            self.active_operations = max(0, self.active_operations - 1)
            self.last_accessed_at = time.time()

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

    def get(self, session_id: str) -> Optional[Session]:
        self._maybe_cleanup()
        now = time.time()
        with self._lock:
            session = self._sessions.get(session_id)
            if session is None:
                return None
            if self._is_expired(session, now):
                del self._sessions[session_id]
                return None
        session.touch()
        return session

    def begin_use(self, session_id: str) -> Optional[Session]:
        """Return a session protected from TTL cleanup during active work."""
        self._maybe_cleanup()
        now = time.time()
        with self._lock:
            session = self._sessions.get(session_id)
            if session is None:
                return None
            if self._is_expired(session, now):
                del self._sessions[session_id]
                return None
            session.begin_use()
            return session

    def end_use(self, session_id: str) -> None:
        """Release an active-operation lease if the session still exists."""
        with self._lock:
            session = self._sessions.get(session_id)
        if session is not None:
            session.end_use()

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
