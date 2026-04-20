#!/bin/bash
# Check for sorry in Lean source code.
# Delegates to check_sorry.py which correctly strips /- ... -/ block
# comments before counting.
#
# Historical note: an earlier bash-only implementation used the regex
#   grep -c '^\s*[^-].*\bsorry\b'
# which matched the word "sorry" inside docstrings that mentioned
# phrases like "sorry-free" or "0 sorry", producing false positives.
# Use the Python version as the source of truth.
#
# Usage:
#   ./scripts/check_sorry.sh           # report all sorry
#   ./scripts/check_sorry.sh --layer1  # exit 1 if Layer 1 has sorry

set -e
cd "$(dirname "$0")/.."

if command -v python3 >/dev/null 2>&1; then
    PY=python3
elif command -v python >/dev/null 2>&1; then
    PY=python
else
    echo "ERROR: python/python3 not found; cannot run sorry audit"
    exit 2
fi

exec "$PY" scripts/check_sorry.py "$@"
