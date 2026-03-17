#!/usr/bin/env python3
"""Run the public TCGA CRC CMS example builder from a source checkout."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from mlheatmap.public_crc_cms import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
