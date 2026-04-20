#!/usr/bin/env python3
"""Tier 7 — Manuscript claim verification (283 checks).

Delegates to verify_manuscript_claims.py which parses the LaTeX
manuscript, supplemental material, and cover letter, then verifies
every numerical claim against authoritative data sources.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

rc = subprocess.call(
    [sys.executable, str(HERE / "verify_manuscript_claims.py")],
    cwd=str(HERE),
)
sys.exit(rc)
