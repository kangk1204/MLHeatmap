"""In-memory session management."""

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
    norm_method: str = ""
    size_factors: Optional[np.ndarray] = None
    gene_names: list = field(default_factory=list)
    sample_names: list = field(default_factory=list)
    groups: dict = field(default_factory=dict)
    excluded_samples: list = field(default_factory=list)
    species: str = "unknown"
    id_type: str = "unknown"
    biomarker_results: Optional[dict] = None
    deg_results: Optional[dict] = None
    heatmap_data: Optional[dict] = None
    heatmap_color_scale: str = "RdBu_r"
    created_at: float = field(default_factory=time.time)


class SessionStore:
    def __init__(self, ttl_hours: int = 4):
        self._sessions: dict[str, Session] = {}
        self._ttl = ttl_hours * 3600
        self._last_cleanup: float = 0.0

    def create(self) -> Session:
        s = Session()
        self._sessions[s.id] = s
        self._maybe_cleanup()
        return s

    def get(self, session_id: str) -> Optional[Session]:
        self._maybe_cleanup()
        session = self._sessions.get(session_id)
        if session and time.time() - session.created_at > self._ttl:
            del self._sessions[session_id]
            return None
        return session

    def _maybe_cleanup(self):
        now = time.time()
        if now - self._last_cleanup > 300:  # every 5 minutes
            self._cleanup()

    def _cleanup(self):
        now = time.time()
        self._last_cleanup = now
        expired = [k for k, v in self._sessions.items() if now - v.created_at > self._ttl]
        for k in expired:
            del self._sessions[k]
