"""
Cross-verification (numerical) of companion manuscript claims.

Uses numpy for matrix operations on density matrices, Pauli operators,
and MLE reconstruction.  Verifies numerical claims from Theorem 1,
Lemma 2, and the experimental predictions.

Requirements:  pip install numpy
Usage:         python numerical_bounds.py

The exact tier covers algebraic claims; this tier covers numerical
predictions requiring state construction and MLE simulation.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root (registry.py)

from test_hessian import test_hessian_properties
from test_mle_basin import test_mle_basin
from test_fidelity_bounds import test_fidelity_bounds
from test_linear_bound import test_linear_bound
from test_independence import test_independence
from test_basis_cover_optimality import test_basis_cover_optimality


_TESTS = [
    ("Hessian properties (Theorem 1, Lemma)", test_hessian_properties),
    ("MLE informative-subspace identifiability (Theorem 1(i))", test_mle_basin),
    ("Fidelity bounds (Lemma 2)", test_fidelity_bounds),
    ("Spectral characterization (Theorem 2)", test_linear_bound),
    ("Independence structure (Lemma 1)", test_independence),
    ("Basis-cover optimality certificate (Sec. compression-advantage)",
     test_basis_cover_optimality),
]

# G1 coverage guard.
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from scripts.discover_tests import assert_coverage as _assert_coverage
_assert_coverage(
    tier_dir=_Path(__file__).resolve().parent,
    registered_modules={fn.__module__ for _, fn in _TESTS},
    pattern="test_*.py",
    ignore={"common", "run_all"},
)


def main():
    print("=" * 68)
    print("CROSS-VERIFICATION: NUMERICAL (numpy)")
    print("companion manuscript: Krawtchouk spectral correspondence")
    print("Tolerances: exact values 1e-10, saturation 1e-8, "
          "eigenstate 1e-12, informative 1e-6")
    print("=" * 68)
    total = 0
    failed_modules = []

    for label, func in _TESTS:
        print(f"\n--- {label} ---")
        try:
            total += func()
        except Exception as exc:
            print(f"  ERROR: {exc}")
            failed_modules.append(label)

    print("\n" + "=" * 68)
    if failed_modules:
        print(f"TOTAL: {total} numerical checks passed, "
              f"{len(failed_modules)} module(s) FAILED:")
        for m in failed_modules:
            print(f"  *** {m}")
    else:
        print(f"TOTAL: {total} numerical checks passed")
    print("=" * 68)
    return total, len(failed_modules)


if __name__ == "__main__":
    total, n_failed = main()
    sys.exit(0 if total > 0 and n_failed == 0 else 1)
