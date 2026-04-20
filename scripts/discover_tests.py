#!/usr/bin/env python3
"""G1 of the CI/CD gap inventory: glob-based test coverage guard.

The existing per-tier `run_all.py` files keep a hand-maintained
`_TESTS = [...]` list that maps (description -> callable).  This guard
*does not replace that list* -- the list still carries the human-readable
description and the optional-per-test entrypoint indirection that
already works.  What it adds is a filesystem check: if a new
`test_*.py` (or configurable pattern) lands in the tier directory but
is not referenced by any entry in the list, we hard-fail loudly
instead of silently dropping the test.

Usage from inside a tier's `run_all.py`:
    from scripts.discover_tests import assert_coverage
    assert_coverage(
        tier_dir=Path(__file__).parent,
        registered_modules={t.__module__ for _, t in _TESTS},
        pattern="test_*.py",
        ignore={"common", "run_all"},
    )

If that import path is inconvenient (tier dirs are not on PYTHONPATH
by default), the function is also importable as a plain script for
ad-hoc auditing:
    python scripts/discover_tests.py tier2-sympy test_*.py
"""
from __future__ import annotations
import argparse
import fnmatch
import sys
from pathlib import Path
from typing import Iterable


def discover_test_modules(tier_dir: Path, pattern: str = "test_*.py",
                          ignore: Iterable[str] = ()) -> list[str]:
    """Return the module names (filename without .py) matching `pattern`
    inside `tier_dir`, excluding any name in `ignore`."""
    ignore_set = set(ignore)
    out: list[str] = []
    for p in sorted(tier_dir.iterdir()):
        if not p.is_file():
            continue
        if not fnmatch.fnmatch(p.name, pattern):
            continue
        mod = p.stem
        if mod in ignore_set or mod.startswith("_"):
            continue
        out.append(mod)
    return out


def assert_coverage(tier_dir: Path, registered_modules: Iterable[str],
                    pattern: str | Iterable[str] = "test_*.py",
                    ignore: Iterable[str] = ()) -> None:
    """Raise SystemExit(2) if a pattern-matching file exists on disk but
    is not in `registered_modules`.  Called at the top of a tier's
    `run_all.py` so that adding a test file without wiring it up is a
    hard error instead of a silent dropout.

    `pattern` may be a single glob (string) or an iterable of globs;
    when multiple, findings from ALL patterns are aggregated before
    exit so the user sees every drift in one run (no iterate-to-find).
    An empty iterable is a coherent "no tests in this tier" state and
    silently passes; if a caller wants the guard to actively assert
    coverage, they must pass at least one pattern.

    Both sides are filtered to each pattern: a tier may register
    multiple conventions (e.g. tier 4 has both `verify_*.py` and
    `test_*.py`), and we only compare within each pattern bucket so
    that registering `test_mle_purity` under `test_*.py` doesn't look
    "stale" when this assertion is called with pattern `verify_*.py`.

    Raises:
        FileNotFoundError: if `tier_dir` does not exist on disk.
        SystemExit(2): if any file matches a pattern but is not
            registered, OR any registered module is absent from disk.
    """
    patterns = [pattern] if isinstance(pattern, str) else list(pattern)
    registered_set = set(registered_modules)
    drift_reports: list[str] = []

    for pat in patterns:
        on_disk = set(discover_test_modules(tier_dir, pat, ignore))
        # Filter registered to names whose `mod + ".py"` would match pattern.
        registered_in_pattern = {
            m for m in registered_set
            if fnmatch.fnmatch(m + ".py", pat)
        }
        missing = sorted(on_disk - registered_in_pattern)
        stale = sorted(registered_in_pattern - on_disk)
        if missing:
            drift_reports.append(
                f"  {len(missing)} file(s) match {pat!r} but are not "
                f"in _TESTS: {missing}"
            )
        if stale:
            drift_reports.append(
                f"  {len(stale)} entry(ies) in _TESTS match {pat!r} but "
                f"reference modules not present on disk: {stale}"
            )

    if drift_reports:
        sys.stderr.write(
            f"\nERROR: test-coverage drift in {tier_dir.name}:\n"
            + "\n".join(drift_reports)
            + "\n  Either add the test to _TESTS or rename the file.\n"
        )
        sys.exit(2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("tier_dir", type=Path)
    ap.add_argument("pattern", nargs="?", default="test_*.py")
    args = ap.parse_args()
    mods = discover_test_modules(args.tier_dir, args.pattern)
    print(f"[discover_tests] {args.tier_dir.name}: {len(mods)} file(s) matching "
          f"{args.pattern!r}:")
    for m in mods:
        print(f"  {m}.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
