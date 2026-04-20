#!/usr/bin/env python3
"""Phase H orchestrator: regenerate every machine-verified artifact.

Runs, in order:
  1. dump_axioms.py                     -> generated/axiom_report.{json,txt}
  2. generate_verification_manifest.py  -> generated/verification_manifest.json
  3. generate_sm_axiom_table.py         -> generated/sm_axiom_table.tex

Exits non-zero if any step fails, so CI can gate on a single command.
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
BASELINE = SCRIPTS.parent / ".baseline" / "verification_manifest.json"
CURRENT = SCRIPTS.parent / "generated" / "verification_manifest.json"

STEPS = [
    ("dump_axioms.py", []),
    ("generate_verification_manifest.py", []),
    ("generate_sm_axiom_table.py", []),
    # Phase I validators (run last; fail-fast on drift).
    ("sync_registry_axioms.py", ["--check"]),
    ("check_generated_links.py", []),
]
# G10 optional tail: manifest-vs-baseline diff.  Only runs when a
# baseline is committed at `.baseline/verification_manifest.json`.
# Paths are passed relative to the working directory (SCRIPTS.parent =
# lean4/) so logs read cleanly on both Windows and POSIX shells.
OPTIONAL_DIFF = ("diff_manifest.py", [
    str(CURRENT.relative_to(SCRIPTS.parent)),
    str(BASELINE.relative_to(SCRIPTS.parent)),
])


def main() -> int:
    steps = list(STEPS)
    if BASELINE.is_file():
        steps.append(OPTIONAL_DIFF)
    for step, extra_args in steps:
        print(f"\n=== running {step} {' '.join(extra_args)} ===")
        proc = subprocess.run(
            [sys.executable, str(SCRIPTS / step), *extra_args],
            cwd=SCRIPTS.parent,
        )
        if proc.returncode != 0:
            print(f"FAIL: {step} exited {proc.returncode}")
            return proc.returncode
    print("\n[generate_all] all artifacts up to date, no drift")
    return 0


if __name__ == "__main__":
    sys.exit(main())
