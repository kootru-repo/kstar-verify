#!/bin/bash
# Entrypoint for the SageMath verification container.  Runs the four
# mathy tiers in the order they compose: SageMath (tier 1) produces
# exact combinatorial certificates, SymPy (tier 2) checks rational
# identities, NumPy (tier 3) checks numerical bounds, and tier 5
# verifies the HS -> F fidelity recovery.  Tiers 4 and 6 need bundled
# hardware data and are run via validate.sh from the repo root.
set -euo pipefail
# Sage's Python buffers stdout by default, so under `docker run` without
# a TTY the entire PASS/FAIL log can be lost if Python exits before the
# buffer flushes.  Force line-buffered stdout so every tier's output
# arrives in the container log as it happens.
export PYTHONUNBUFFERED=1

echo "=== Tier 1: SageMath exact arithmetic (44 tests) ==="
cd /verify/tier1-sage
sage run_all.sage --export
echo

echo "=== Tier 2: SymPy exact rational arithmetic ==="
cd /verify/tier2-sympy
sage -python run_all.py
echo

echo "=== Tier 3: NumPy numerical bounds ==="
cd /verify/tier3-numpy
sage -python run_all.py
echo

echo "=== Tier 5: Fidelity HS->F recovery (10 checks) ==="
cd /verify/tier5-fidelity
sage -python test_fidelity_recovery.py
echo

echo "Tiers 1, 2, 3, 5 passed (tiers 4, 6 require hardware data; use validate.sh)."
