#!/usr/bin/env python3
"""
Run all independent verification scripts and report summary.

Runs scripts in parallel (independent of each other) and collects results.
"""
import subprocess
import sys
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

SCRIPTS = [
    ("SM Prop 1 & Lemma 6 (SymPy exact)", "verify_propositions.py"),
    ("K* operator set & compression", "verify_operator_set.py"),
    ("Hardware fidelities (cvxpy + robust_mle)", "verify_hardware_fidelities.py"),
    ("GHZ reanalysis (DFE + subspace)", "verify_ghz_reanalysis.py"),
    ("Allocation fractions & D(S-AR)", "verify_allocation.py"),
    ("SOTA comparison, statistics & decomposition", "verify_sota_and_stats.py"),
    ("RDM fidelities (partial trace + lstsq)", "verify_rdm_fidelities.py"),
    ("Qutrit extension (Gell-Mann + q-ary Krawtchouk)", "verify_qutrit.py"),
    ("Fisher information, frame theory & ablation", "verify_fisher_frame.py"),
    ("Figure data consistency (gen_figures vs JSON)", "verify_figures.py"),
    ("Extended SM tables (qudit, ordering, scaling)", "verify_sm_extended.py"),
    ("MLE purity (Thm 1(iii) near-purity condition)", "test_mle_purity.py"),
]

# Slow simulation-only scripts (run with --full flag)
SLOW_SCRIPTS = [
    ("Remaining claims (ablation, hedge, convergence)", "verify_remaining_claims.py"),
]

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent  # kstar-verify root (registry.py lives here)
TIMEOUT = 1200  # seconds per script (percentile + remaining need >600s)

# G1 coverage guard: hard-fail if a verify_*.py or test_*.py file exists
# on disk but is not in SCRIPTS + SLOW_SCRIPTS.  Tier 4 uses a filename-
# based list (not importlib), so we convert both sides to module names.
# Pass BOTH patterns in one call so every drift is reported at once
# instead of requiring a second run to see the test_*.py bucket.
sys.path.insert(0, str(REPO_ROOT))
from scripts.discover_tests import assert_coverage as _assert_coverage
_registered = {Path(script).stem for _, script in SCRIPTS + SLOW_SCRIPTS}
_assert_coverage(
    tier_dir=HERE,
    registered_modules=_registered,
    pattern=("verify_*.py", "test_*.py"),
    ignore={"common", "run_all", "core", "robust_mle"},
)

# Auto-discover bundled data/ if KSTAR_DATA_DIR is not set
if "KSTAR_DATA_DIR" not in os.environ:
    default_data = REPO_ROOT / "data"
    if default_data.is_dir():
        os.environ["KSTAR_DATA_DIR"] = str(default_data)

# Ensure child processes can import registry.py from repo root
_pypath = os.environ.get("PYTHONPATH", "")
os.environ["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + _pypath if _pypath else "")


def run_one(name, script):
    """Run a single verification script, return (name, script, passed, rc, output, error)."""
    path = HERE / script
    if not path.exists():
        return (name, script, False, -1,
                f"  Script not found: {path}\n", "FileNotFoundError")
    try:
        proc = subprocess.run(
            [sys.executable, str(path)],
            cwd=str(HERE),
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
        )
        return (name, script, proc.returncode == 0, proc.returncode,
                proc.stdout, proc.stderr)
    except subprocess.TimeoutExpired:
        return (name, script, False, -2,
                f"  TIMEOUT after {TIMEOUT}s\n", "TimeoutExpired")
    except OSError as e:
        return (name, script, False, -3,
                f"  OS error: {e}\n", str(e))


def main():
    full = "--full" in sys.argv
    parallel = "--sequential" not in sys.argv
    t0 = time.time()

    scripts = SCRIPTS + (SLOW_SCRIPTS if full else [])

    print("=" * 70)
    print("  INDEPENDENT VERIFICATION SUITE")
    print(f"  Mode: {'parallel' if parallel else 'sequential'}"
          f"{', full (incl. slow sims)' if full else ''}")
    print(f"  Scripts: {len(scripts)}" + (f" (pass --full for {len(SLOW_SCRIPTS)} more)" if not full else ""))
    print("=" * 70)

    results = []

    if parallel:
        # Launch all scripts concurrently
        futures = {}
        with ProcessPoolExecutor(max_workers=len(SCRIPTS)) as pool:
            for name, script in scripts:
                fut = pool.submit(run_one, name, script)
                futures[fut] = (name, script)

            for fut in as_completed(futures):
                results.append(fut.result())

        # Sort back to original order
        order = {s: i for i, (_, s) in enumerate(SCRIPTS)}
        results.sort(key=lambda r: order.get(r[1], 99))
    else:
        for name, script in scripts:
            results.append(run_one(name, script))

    # Print captured output per script
    for name, script, passed, rc, stdout, stderr in results:
        print(f"\n{'='*70}")
        print(f"  Script: {script}  ({'PASS' if passed else 'FAIL'})", flush=True)
        print(f"{'='*70}")
        if stdout:
            print(stdout, end="", flush=True)
        if not passed and stderr:
            # On failure, surface the FULL stderr (not just the tail) so
            # CI logs contain the assertion message / traceback that
            # caused exit != 0.  Truncating to 5 lines historically hid
            # the cause when the last lines were progress prints.
            print(f"  --- STDERR ({script}) ---", flush=True)
            for line in stderr.strip().split("\n"):
                print(f"  STDERR: {line}", flush=True)
            print(f"  --- end STDERR ({script}) ---", flush=True)

    # Summary
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    all_pass = True
    for name, script, passed, rc, _, _ in results:
        if rc == -1:
            tag = "FAIL (not found)"
        elif rc == -2:
            tag = "FAIL (timeout)"
        elif rc == -3:
            tag = "FAIL (OS error)"
        else:
            tag = "PASS" if passed else f"FAIL (exit {rc})"
        print(f"  [{tag}] {name}")
        if not passed:
            all_pass = False

    print(f"\n  Elapsed: {elapsed:.1f}s")
    if all_pass:
        print("  ALL VERIFICATION SUITES PASSED")
    else:
        print("  *** SOME SUITES FAILED ***")
    print("=" * 70)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
