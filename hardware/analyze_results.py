"""Offline analyzer for live-QPU W-state K* JSON output.

Pointed at a results JSON produced by run_w_state.py (or at the
shipped reference data/w_repeat_results.json), this prints the
per-run table, mean +/- std, a bootstrap 95 % CI on dF, and a
bound-check against hardware/expected_bounds.json.

Usage:
  python hardware/analyze_results.py <results.json>
  python hardware/analyze_results.py --reference
  python hardware/analyze_results.py --all hardware/results/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
REFERENCE_JSON = REPO_ROOT / "data" / "w_repeat_results.json"
BOUNDS_JSON = THIS_DIR / "expected_bounds.json"


def _print_table(runs: list[dict]) -> None:
    print(f"  {'run':>4s}  {'seed':>6s}  {'F(K*)':>8s}  {'F(rand)':>8s}  {'dF':>8s}")
    print(f"  {'-' * 4}  {'-' * 6}  {'-' * 8}  {'-' * 8}  {'-' * 8}")
    for r in runs:
        print(f"  {r['run']:>4d}  {r['seed']:>6d}  "
              f"{r['f_kstar']:>8.4f}  {r['f_rand']:>8.4f}  {r['delta_f']:>+8.4f}")


def analyze(path: Path) -> int:
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    backend = data.get("backend", "(unknown)")
    runs = data.get("runs", [])
    if not runs:
        print(f"  No runs found in {path}")
        return 2

    print(f"\n  File:    {path}")
    print(f"  Backend: {backend}")
    print(f"  Runs:    {len(runs)} x {data.get('n_shots', '?')} shots")
    print()
    _print_table(runs)

    f_kstars = np.array([r["f_kstar"] for r in runs])
    f_rands = np.array([r["f_rand"] for r in runs])
    deltas = np.array([r["delta_f"] for r in runs])
    ddof = 1 if len(runs) > 1 else 0

    print(f"\n  F(K*)   mean = {f_kstars.mean():.4f}   std = {f_kstars.std(ddof=ddof):.4f}")
    print(f"  F(rand) mean = {f_rands.mean():.4f}   std = {f_rands.std(ddof=ddof):.4f}")
    print(f"  dF      mean = {deltas.mean():+.4f}   std = {deltas.std(ddof=ddof):.4f}")

    if len(runs) >= 2:
        rng = np.random.default_rng(0)
        boot = rng.choice(deltas, size=(10_000, len(deltas)), replace=True).mean(axis=1)
        print(f"  dF 95% CI = [{np.percentile(boot, 2.5):+.4f}, "
              f"{np.percentile(boot, 97.5):+.4f}]   "
              f"(bootstrap, 10000 resamples)")

    return _check_bounds(f_kstars, deltas, len(runs))


def _check_bounds(f_kstars: np.ndarray, deltas: np.ndarray, n_runs: int) -> int:
    if not BOUNDS_JSON.is_file():
        print("\n  No expected_bounds.json; skipping bound check.")
        return 0
    with BOUNDS_JSON.open(encoding="utf-8") as f:
        b = json.load(f)

    checks = [
        ("F(K*) mean", f_kstars.mean() >= b.get("f_kstar_min", 0.0),
         f_kstars.mean(), b.get("f_kstar_min", 0.0)),
        ("dF mean",     deltas.mean() >= b.get("delta_f_min", 0.0),
         deltas.mean(), b.get("delta_f_min", 0.0)),
        ("runs with dF > 0", (deltas > 0).sum() >= b.get("n_positive_runs_min", 0),
         int((deltas > 0).sum()), b.get("n_positive_runs_min", 0)),
    ]
    print("\n  Bound check (hardware/expected_bounds.json):")
    all_ok = True
    for name, ok, got, need in checks:
        tag = "PASS" if ok else "FAIL"
        all_ok &= ok
        print(f"    [{tag}] {name}: got {got}, need >= {need}")
    return 0 if all_ok else 1


def analyze_all(dir_path: Path) -> int:
    jsons = sorted(dir_path.glob("w_repeat_results_*.json"))
    if not jsons:
        print(f"  No w_repeat_results_*.json under {dir_path}")
        return 2
    all_deltas: list[float] = []
    overall_ok = True
    for p in jsons:
        rc = analyze(p)
        overall_ok &= (rc == 0)
        with p.open(encoding="utf-8") as f:
            data = json.load(f)
        all_deltas.extend(r["delta_f"] for r in data.get("runs", []))
    if all_deltas:
        arr = np.array(all_deltas)
        print(f"\n  COMBINED across {len(jsons)} files "
              f"({len(arr)} runs total): dF mean = {arr.mean():+.4f}   "
              f"std = {arr.std(ddof=1):.4f}")
    return 0 if overall_ok else 1


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group()
    g.add_argument("path", nargs="?", help="Results JSON to analyze")
    g.add_argument("--reference", action="store_true",
                   help="Analyze the shipped reference at data/w_repeat_results.json")
    g.add_argument("--all", dest="all_dir",
                   help="Analyze every w_repeat_results_*.json under a directory")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.reference:
        return analyze(REFERENCE_JSON)
    if args.all_dir:
        return analyze_all(Path(args.all_dir))
    if not args.path:
        print("ERROR: provide a path, or pass --reference / --all DIR")
        return 2
    return analyze(Path(args.path))


if __name__ == "__main__":
    raise SystemExit(main())
