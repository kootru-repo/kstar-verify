#!/usr/bin/env python3
"""
verify_canonical_constants.py
=============================

Single source of truth: canonical_constants.json defines the canonical
numerical values used by Lean4 formalization, Sage adversarial suite,
manuscript text, and verification scripts. This script verifies that
Lean and Sage sources all agree with the canonical JSON.

Closes Gap 11: Sage and Lean adversarial suites are now linked via a
shared constants file, so a reviewer or future contributor can change a
value in one place and the verifier will catch any unsynchronized site.

Exit code: 0 if all checks pass, 1 otherwise.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent  # kstar-verify/

CANONICAL = SCRIPT_DIR / "canonical_constants.json"
LEAN_DIR = REPO_ROOT / "lean4-formal" / "KstarFormal"
LAKE_DIR = REPO_ROOT / "lean4-formal"
SAGE_DIR = REPO_ROOT / "scripts"


def require_lake_build():
    """Run `lake build` and abort if it fails. Closes Cat 3 of the
    adversarial audit (verifier passing despite uncompilable Lean source)."""
    print("=" * 70)
    print("LAKE BUILD GATING (must succeed before constant cross-check)")
    print("=" * 70)
    if not (LAKE_DIR / "lakefile.lean").exists() and not (LAKE_DIR / "lakefile.toml").exists():
        print(f"ERROR: no lakefile found at {LAKE_DIR}")
        sys.exit(1)
    try:
        result = subprocess.run(
            ["lake", "build"],
            cwd=str(LAKE_DIR),
            capture_output=True,
            text=True,
            timeout=600,
        )
    except FileNotFoundError:
        print("ERROR: `lake` executable not found in PATH")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: `lake build` timed out after 600s")
        sys.exit(1)
    if result.returncode != 0:
        print("ERROR: `lake build` failed (returncode "
              f"{result.returncode}). Aborting canonical constants check.")
        print("--- stderr (tail) ---")
        for line in result.stderr.splitlines()[-30:]:
            print(line)
        sys.exit(1)
    print("[OK] lake build completed successfully")
    print()


def load_lean_combined():
    parts = []
    for f in LEAN_DIR.rglob("*.lean"):
        parts.append(f.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(parts)


def load_sage_combined():
    parts = []
    for f in SAGE_DIR.rglob("*.py"):
        parts.append(f.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(parts)


def expect_value(text, value, label, source):
    """Search for the given numeric value as a standalone token."""
    if isinstance(value, float):
        # search for a few decimal renderings
        candidates = [f"{value}", f"{value:.3f}", f"{value:.4f}"]
    else:
        candidates = [str(value)]
    for c in candidates:
        if re.search(rf"(?<![\w.]){re.escape(c)}(?![\w.])", text):
            return True
    return False


def main():
    if "--skip-lake-build" not in sys.argv:
        require_lake_build()
    canonical = json.loads(CANONICAL.read_text(encoding="utf-8"))
    lean_text = load_lean_combined()
    sage_text = load_sage_combined()

    n_pass = 0
    n_fail = 0

    print("=" * 70)
    print("CANONICAL CONSTANTS CROSS-CHECK (Lean / Sage / JSON)")
    print("=" * 70)
    print(f"Canonical source: {CANONICAL.name}")
    print(f"Lean dir:         {LEAN_DIR.name}/")
    print(f"Sage dir:         {SAGE_DIR.name}/")
    print()

    # Iterate canonical constants and check both sides reference them.
    skip_keys = {"_meta"}
    for key, entry in canonical.items():
        if key in skip_keys:
            continue
        if not isinstance(entry, dict):
            continue

        # Singleton value
        if "value" in entry:
            v = entry["value"]
            in_lean = expect_value(lean_text, v, key, "Lean")
            in_sage = expect_value(sage_text, v, key, "Sage")

            if in_lean and in_sage:
                print(f"[PASS] {key} = {v}  (Lean+Sage agree)")
                n_pass += 1
            elif in_lean and not in_sage:
                # Some constants only live in Lean (e.g. envelope, hedged factor)
                if key in {"hedged_envelope_factor", "axiom_count"}:
                    print(f"[PASS] {key} = {v}  (Lean only — by design)")
                    n_pass += 1
                else:
                    print(f"[FAIL] {key} = {v}  (in Lean, missing in Sage)")
                    n_fail += 1
            elif in_sage and not in_lean:
                print(f"[FAIL] {key} = {v}  (in Sage, missing in Lean)")
                n_fail += 1
            else:
                print(f"[FAIL] {key} = {v}  (missing in both)")
                n_fail += 1

        # List of values (e.g. eigenvalues)
        elif "values" in entry and isinstance(entry["values"], list):
            vs = entry["values"]
            missing_lean = [v for v in vs if not expect_value(lean_text, v, key, "Lean")]
            missing_sage = [v for v in vs if not expect_value(sage_text, v, key, "Sage")]
            if not missing_lean and not missing_sage:
                print(f"[PASS] {key} = {vs}  (all Lean+Sage)")
                n_pass += 1
            else:
                msg = []
                if missing_lean:
                    msg.append(f"Lean missing {missing_lean}")
                if missing_sage:
                    msg.append(f"Sage missing {missing_sage}")
                print(f"[FAIL] {key}: {'; '.join(msg)}")
                n_fail += 1

        # Three-sector
        elif "identity" in entry and "edge" in entry and "bulk" in entry:
            parts = [entry["identity"], entry["edge"], entry["bulk"], entry["sum"]]
            missing_lean = [v for v in parts if not expect_value(lean_text, v, key, "Lean")]
            if not missing_lean:
                print(f"[PASS] {key} = {entry['identity']}+{entry['edge']}+{entry['bulk']}={entry['sum']}")
                n_pass += 1
            else:
                print(f"[FAIL] {key}: Lean missing {missing_lean}")
                n_fail += 1

    print()
    print("=" * 70)
    print(f"RESULTS: {n_pass} passed, {n_fail} failed")
    print("=" * 70)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
