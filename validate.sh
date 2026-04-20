#!/usr/bin/env bash
# K* Verification Suite — one-command validation
# Usage:
#   ./validate.sh              # Full: Lean4 + SageMath (Docker) + Python tiers
#   ./validate.sh --quick      # Python only (no Docker, no Lean)
#   ./validate.sh --no-lean    # Skip Lean4 build (Docker + Python)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
# Docker on Windows+MSYS needs backslash-free absolute paths; cygpath -m
# converts MSYS-style /c/foo to mixed C:/foo. Only used for the Tier 1
# docker mount; native bash paths are kept everywhere else so the Python
# tier runners on Windows still receive proper Windows-style paths.
if command -v cygpath >/dev/null 2>&1; then
    ROOT_DOCKER="$(cygpath -m "$ROOT")"
else
    ROOT_DOCKER="$ROOT"
fi
QUICK=false
NO_LEAN=false
FAILED=0

# Detect Python: prefer python3 if it has numpy, else fall back to python
PYTHON=python3
if ! "$PYTHON" -c "import numpy" 2>/dev/null; then
    PYTHON=python
    if ! "$PYTHON" -c "import numpy" 2>/dev/null; then
        echo "ERROR: neither python3 nor python has numpy installed."
        echo "Run: pip install -r $ROOT/requirements.txt"
        exit 1
    fi
fi

for arg in "$@"; do
    case "$arg" in
        --quick)   QUICK=true; NO_LEAN=true ;;
        --no-lean) NO_LEAN=true ;;
    esac
done

run_tier() {
    local name="$1"; shift
    echo ""
    echo "========================================"
    echo "  $name"
    echo "========================================"
    if "$@"; then
        echo "  PASS"
    else
        echo "  FAIL"
        FAILED=$((FAILED + 1))
    fi
}

# --- Lean4 (Tier 0) ---
if [ "$NO_LEAN" = false ]; then
    run_tier "Tier 0: Lean4 formal proofs (0 sorry)" \
        bash -c "cd '$ROOT/lean4' && lake build"
fi

# --- SageMath via Docker (Tier 1) ---
# Mount the full kstar-verify root so the Sage suite can reach the
# Lean source (lean4/KstarFormal/Combinatorics/GhzNonCoverage.lean) for
# the registry-fact witnesses that parse the Lean constant. Working
# directory is the tier1-sage subtree.
if [ "$QUICK" = false ]; then
    run_tier "Tier 1: SageMath exact arithmetic (44 tests)" \
        env MSYS_NO_PATHCONV=1 docker run --rm \
            -v "${ROOT_DOCKER}:/work" \
            -w /work/tier1-sage \
            sagemath/sagemath:latest \
            sage run_all.sage --export
fi

# --- Python tiers (2, 3, 4, 5, 6) ---
run_tier "Tier 2: SymPy exact rational arithmetic" \
    "$PYTHON" "$ROOT/tier2-sympy/run_all.py"

run_tier "Tier 3: NumPy numerical bounds" \
    "$PYTHON" "$ROOT/tier3-numpy/run_all.py"

# Auto-discover bundled data if KSTAR_DATA_DIR is not set
if [ -z "${KSTAR_DATA_DIR:-}" ] && [ -d "$ROOT/data" ]; then
    export KSTAR_DATA_DIR="$ROOT/data"
    echo "  (auto-discovered bundled data at $ROOT/data)"
fi

if [ -n "${KSTAR_DATA_DIR:-}" ]; then
    run_tier "Tier 4: Independent verification (data-consistency)" \
        "$PYTHON" "$ROOT/tier4-independent/run_all.py"
else
    echo ""
    echo "  Tier 4: SKIPPED (data/ not found under $ROOT; re-clone to restore)"
fi

run_tier "Tier 5: Fidelity HS->F recovery (10 checks)" \
    "$PYTHON" "$ROOT/tier5-fidelity/test_fidelity_recovery.py"

if [ -n "${KSTAR_DATA_DIR:-}" ]; then
    run_tier "Tier 6: End-to-end implication chains" \
        "$PYTHON" "$ROOT/tier6-chains/test_implication_chain.py"
else
    echo ""
    echo "  Tier 6: SKIPPED (data/ not found under $ROOT; re-clone to restore)"
fi

# --- Summary ---
echo ""
echo "========================================"
if [ "$FAILED" -eq 0 ]; then
    echo "  ALL TIERS PASSED"
else
    echo "  $FAILED TIER(S) FAILED"
    exit 1
fi
