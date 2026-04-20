"""
Cross-verification (exact arithmetic) of companion manuscript proof chains.

Uses SymPy for exact integer/rational arithmetic.  Every assertion maps
to a specific theorem, proposition, or claim in the manuscript.

Requirements:  pip install sympy
Usage:         python exact_arithmetic.py
               python exact_arithmetic.py --sequential

All functions use SymPy exact arithmetic (ZZ/QQ) — no floating point
in any algebraic claim.
"""

import sys
import os
import io
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root (registry.py)

import test_krawtchouk
import test_gram_krawtchouk
import test_basin_separation
import test_purity_bound
import test_universality
import test_weight_saturation
import test_asymptotic
import test_audit_remediation
import test_gram_full_rank
import test_bose_mesner_iff

_TESTS = [
    ("Krawtchouk polynomials (Prop 1, SM S5)", test_krawtchouk.test_krawtchouk),
    ("Gram-Krawtchouk diagonalisation (Prop 1)", test_gram_krawtchouk.test_gram_krawtchouk_diagonalisation),
    ("Basin Separation Theorem (Theorem 1)", test_basin_separation.test_basin_separation),
    ("Purity Bound (Lemma 2)", test_purity_bound.test_purity_bound),
    ("Universality across Hamming schemes (Lemma 5)", test_universality.test_universality),
    ("Weight saturation properties (Theorem 1 precondition)", test_weight_saturation.test_weight_saturation),
    ("Asymptotic scaling (Theorem 3)", test_asymptotic.test_asymptotic),
    ("Audit remediation (Findings 5, 7, 11)", test_audit_remediation.main),
    ("Gram matrix full rank G(n) (Prop 1)", test_gram_full_rank.test_gram_full_rank),
    ("Bose-Mesner iff characterization (Lemma 5)", test_bose_mesner_iff.test_bose_mesner_iff),
]

# G1 coverage guard: hard-fail if a test_*.py file exists on disk but is
# not referenced in _TESTS (or vice versa).  Prevents silent dropouts.
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from scripts.discover_tests import assert_coverage as _assert_coverage
_assert_coverage(
    tier_dir=_Path(__file__).resolve().parent,
    registered_modules={fn.__module__ for _, fn in _TESTS},
    pattern="test_*.py",
    ignore={"common", "run_all"},
)


def _run_test_worker(qualified_name):
    """Run a test function in a worker process, capturing stdout.

    qualified_name is 'module_name.func_name' (e.g. 'test_krawtchouk.test_krawtchouk').
    """
    import importlib
    _dir = os.path.dirname(os.path.abspath(__file__))
    _root = os.path.dirname(_dir)
    if _dir not in sys.path:
        sys.path.insert(0, _dir)
    if _root not in sys.path:
        sys.path.insert(0, _root)
    mod_name, func_name = qualified_name.rsplit(".", 1)
    mod = importlib.import_module(mod_name)
    func = getattr(mod, func_name)

    buf = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    try:
        count = func()
    finally:
        sys.stdout = old_stdout
    return count, buf.getvalue()


def main():
    ncpu = os.cpu_count() or 4
    nworkers = min(ncpu, len(_TESTS))

    print("=" * 68)
    print("CROSS-VERIFICATION: EXACT ARITHMETIC (SymPy)")
    print("companion manuscript: Krawtchouk spectral correspondence")
    print(f"Parallel execution: {nworkers} workers on {ncpu} cores")
    print("=" * 68)

    results = {}
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        futures = {
            executor.submit(
                _run_test_worker,
                f"{func.__module__}.{func.__name__}"
            ): f"{func.__module__}.{func.__name__}"
            for _, func in _TESTS
        }
        for future in as_completed(futures):
            func_name = futures[future]
            try:
                count, output = future.result()
                results[func_name] = (count, output, False)
            except Exception as exc:
                print(f"\n  FAILED: {func_name}: {exc}")
                results[func_name] = (0, f"  ERROR: {exc}\n", True)

    total = 0
    failed_modules = []
    for label, func in _TESTS:
        entry = results[f"{func.__module__}.{func.__name__}"]
        count, output = entry[0], entry[1]
        is_error = len(entry) > 2 and entry[2]
        print(f"\n--- {label} ---")
        print(output, end="")
        total += count
        if is_error:
            failed_modules.append(label)

    print("\n" + "=" * 68)
    if failed_modules:
        print(f"TOTAL: {total} exact-arithmetic checks passed, "
              f"{len(failed_modules)} module(s) FAILED:")
        for m in failed_modules:
            print(f"  *** {m}")
    else:
        print(f"TOTAL: {total} exact-arithmetic checks passed")
    print("=" * 68)
    return total, len(failed_modules)


def main_sequential():
    """Sequential fallback."""
    print("=" * 68)
    print("CROSS-VERIFICATION: EXACT ARITHMETIC (SymPy)")
    print("companion manuscript: Krawtchouk spectral correspondence")
    print("Sequential execution")
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
        print(f"TOTAL: {total} exact-arithmetic checks passed, "
              f"{len(failed_modules)} module(s) FAILED:")
        for m in failed_modules:
            print(f"  *** {m}")
    else:
        print(f"TOTAL: {total} exact-arithmetic checks passed")
    print("=" * 68)
    return total, len(failed_modules)


if __name__ == "__main__":
    if "--sequential" in sys.argv:
        total, n_failed = main_sequential()
    else:
        total, n_failed = main()
    sys.exit(0 if total > 0 and n_failed == 0 else 1)
