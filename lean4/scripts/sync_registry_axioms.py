#!/usr/bin/env python3
"""Phase I: sync `lean4_axioms` and `lean4_universal_axioms` fields
in proofs_registry.yaml from the machine-verified axiom_report.json.

Only the NAMED propositional axioms are written into the registry
(Lean core, native_decide certificates, and opaque declarations are
noise for the registry's trust-base view).

Idempotent: exits 0 with "already in sync" if nothing changes.
Writes a backup file `proofs_registry.yaml.bak` before modifying.

Usage:
  python scripts/sync_registry_axioms.py [--dry-run] [--check]

  --dry-run:  print the diff without writing
  --check:    exit 1 if registry is out of sync (for CI)
"""
from __future__ import annotations
import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

try:
    import yaml
except ModuleNotFoundError:
    sys.stderr.write(
        "ERROR: PyYAML required.  Install: pip install pyyaml>=6.0\n"
    )
    sys.exit(2)

LEAN4_ROOT = Path(__file__).resolve().parent.parent
KSTAR_ROOT = LEAN4_ROOT.parent
REGISTRY = KSTAR_ROOT / "proofs_registry.yaml"
AXIOM_REPORT = LEAN4_ROOT / "generated" / "axiom_report.json"

NAMED_PROPOSITIONAL_AXIOMS = {
    "witness_states_are_valid",
    "fidelity_orthogonal_zero",
    "fidelity_witness_traceproj_decomp",
    "fidelity_nonneg",
    "weyl_eigenvalue_bound",
}


def named_only(axioms: list[str]) -> list[str]:
    return sorted({a for a in axioms if a in NAMED_PROPOSITIONAL_AXIOMS})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    if not AXIOM_REPORT.is_file():
        print(f"ERROR: {AXIOM_REPORT} missing. Run dump_axioms.py first.")
        return 2

    axiom_report = json.loads(AXIOM_REPORT.read_text(encoding="utf-8"))
    raw = REGISTRY.read_text(encoding="utf-8")
    data = yaml.safe_load(raw)

    changes: list[str] = []
    for stmt in data.get("statements", []):
        cid = stmt.get("id", "<unknown>")
        for reg_key, lean_key in (
            ("lean4_axioms", "lean4_theorem"),
            ("lean4_universal_axioms", "lean4_theorem_universal"),
        ):
            tname = stmt.get(lean_key)
            if not tname or tname not in axiom_report:
                continue
            derived = named_only(axiom_report[tname])
            current = stmt.get(reg_key)
            if current is None:
                # Field absent: only add for universal entries that have
                # a universal theorem registered.
                if reg_key == "lean4_universal_axioms":
                    stmt[reg_key] = derived
                    changes.append(f"+ {cid}.{reg_key} = {derived}")
                else:
                    if derived != []:
                        stmt[reg_key] = derived
                        changes.append(f"+ {cid}.{reg_key} = {derived}")
                    else:
                        stmt[reg_key] = []
                        changes.append(f"+ {cid}.{reg_key} = []")
                continue
            if not isinstance(current, list):
                # Registry field is malformed (e.g. a string).  Treat as
                # out-of-sync; overwrite with the canonical list form.
                changes.append(
                    f"~ {cid}.{reg_key}: non-list {current!r} -> {derived}"
                )
                stmt[reg_key] = derived
                continue
            if sorted(current) != derived:
                changes.append(
                    f"~ {cid}.{reg_key}: {sorted(current)} -> {derived}"
                )
                stmt[reg_key] = derived

    if not changes:
        print("[sync_registry] already in sync (no changes)")
        return 0

    if args.check:
        print(f"[sync_registry] OUT OF SYNC: {len(changes)} change(s)")
        for c in changes[:20]:
            print(f"  {c}")
        if len(changes) > 20:
            print(f"  ... and {len(changes)-20} more")
        return 1

    print(f"[sync_registry] {len(changes)} change(s):")
    for c in changes:
        print(f"  {c}")

    if args.dry_run:
        return 0

    backup = REGISTRY.with_suffix(".yaml.bak")
    shutil.copy2(REGISTRY, backup)
    # Atomic write: temp file in the same directory, then os.replace.
    # If Ctrl-C fires mid-write, the registry is untouched (the temp
    # file is orphaned and can be ignored; it starts with a dot and
    # ends with .tmp, so it's easy to spot and clean up).
    payload = yaml.safe_dump(data, sort_keys=False, width=120)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(REGISTRY.parent), prefix=f".{REGISTRY.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
            f.write(payload)
        os.replace(tmp_name, REGISTRY)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise
    print(f"[sync_registry] wrote {REGISTRY}  (backup: {backup})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
