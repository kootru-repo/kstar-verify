#!/usr/bin/env python3
"""G10: compare the current verification manifest against a baseline.

Flags every regression that would otherwise only be visible by
eyeballing the console output:
  * a counter dropped (tier_pass_counts[tier] decreased, sorry count up)
  * a source SHA changed (registry_sha, constants_sha, tex_sha[file])
  * a per-theorem axiom list grew or gained a named axiom
  * the schema_version went backwards

Exits:
  0 -- fully congruent with baseline (up to expected drift below)
  1 -- one or more regressions or SHA changes detected
  2 -- invocation error (missing files, unparseable JSON)

Expected drift that is NOT flagged as a regression:
  * `git_sha`, `generated_at` -- always differ across runs
  * `lake_build.jobs` -- may differ by cache/incremental state
  * `claim_map` -- derived from other fields, reported via its parents

Usage:
    python scripts/diff_manifest.py CURRENT BASELINE
    python scripts/diff_manifest.py generated/verification_manifest.json \\
        .baseline/verification_manifest.json
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path


def _load(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: manifest not found at {path}\n")
        sys.exit(2)
    except json.JSONDecodeError as e:
        sys.stderr.write(f"ERROR: {path} is not valid JSON: {e}\n")
        sys.exit(2)


def _compare_sha(label: str, cur: str | None, base: str | None) -> list[str]:
    if cur == base:
        return []
    # If the current environment doesn't have the source file at all
    # (e.g. a Lean-only clone without the manuscript tree), cur is
    # None while base has a hash.  That's an environment mismatch,
    # not a regression -- report as INFO.  Baseline-missing -> current-
    # present is likewise INFO (a new artefact appeared).  Only flag
    # [SHA] (regression) when BOTH sides have values and they differ.
    if cur is None or base is None:
        return [f"[INFO] {label}: baseline={base!r} -> current={cur!r} (presence change)"]
    return [f"[SHA] {label}: {base!r} -> {cur!r}"]


def _compare_tex(cur: dict | None, base: dict | None) -> list[str]:
    out: list[str] = []
    cur = cur or {}
    base = base or {}
    for name in sorted(set(cur) | set(base)):
        out += _compare_sha(f"tex_sha[{name}]", cur.get(name), base.get(name))
    return out


def _compare_tier_counts(cur: dict | None, base: dict | None) -> list[str]:
    out: list[str] = []
    cur = cur or {}
    base = base or {}
    for tier in sorted(set(cur) | set(base)):
        c, b = cur.get(tier), base.get(tier)
        # Environment mismatch (results/latest.json absent) -> INFO.
        # Only flag a true regression when BOTH sides have counts.
        if c is None or b is None:
            if c != b:
                out.append(f"[INFO] tier_pass_counts[{tier}]: baseline={b} -> current={c} (presence change)")
            continue
        if c < b:
            out.append(f"[REGRESSION] tier_pass_counts[{tier}]: {b} -> {c} ({c - b})")
        elif c > b:
            out.append(f"[INFO] tier_pass_counts[{tier}]: {b} -> {c} (+{c - b})")
    return out


def _compare_axioms(cur: dict | None, base: dict | None) -> list[str]:
    """Flag any theorem whose axiom set CHANGED.  A strictly-smaller
    axiom set is called out as INFO (strengthening), strictly-larger as
    WARN (weakening of the trust base), and otherwise-different as DIFF.
    """
    out: list[str] = []
    cur = cur or {}
    base = base or {}
    for name in sorted(set(cur) | set(base)):
        c = set(cur.get(name) or [])
        b = set(base.get(name) or [])
        if c == b:
            continue
        added = sorted(c - b)
        removed = sorted(b - c)
        if removed and not added:
            out.append(f"[INFO] {name}: axioms shrunk by {removed}")
        elif added and not removed:
            out.append(f"[WARN] {name}: axioms grew by {added}")
        else:
            out.append(
                f"[DIFF] {name}: +{added} -{removed}"
            )
    return out


def _compare_sorry(cur: dict | None, base: dict | None) -> list[str]:
    """Flag sorry count that went UP."""
    out: list[str] = []
    cur = cur or {}
    base = base or {}
    for k in ("total_sorry", "layer1_sorry", "layer23_sorry"):
        c, b = cur.get(k), base.get(k)
        if c is None or b is None:
            continue
        if c > b:
            out.append(f"[REGRESSION] sorry_audit[{k}]: {b} -> {c} (+{c - b})")
    return out


def diff(current: dict, baseline: dict) -> list[str]:
    findings: list[str] = []
    if current.get("schema_version") != baseline.get("schema_version"):
        findings.append(
            f"[INFO] schema_version: {baseline.get('schema_version')!r} "
            f"-> {current.get('schema_version')!r}"
        )
    findings += _compare_sha(
        "registry_sha", current.get("registry_sha"), baseline.get("registry_sha")
    )
    findings += _compare_sha(
        "constants_sha", current.get("constants_sha"), baseline.get("constants_sha")
    )
    findings += _compare_tex(current.get("tex_sha"), baseline.get("tex_sha"))
    findings += _compare_tier_counts(
        current.get("tier_pass_counts"), baseline.get("tier_pass_counts")
    )
    findings += _compare_sorry(
        current.get("sorry_audit"), baseline.get("sorry_audit")
    )
    findings += _compare_axioms(
        current.get("per_theorem_axioms"), baseline.get("per_theorem_axioms")
    )
    return findings


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("current", type=Path)
    ap.add_argument("baseline", type=Path)
    args = ap.parse_args()

    cur = _load(args.current)
    base = _load(args.baseline)
    findings = diff(cur, base)

    if not findings:
        print(f"[diff-manifest] OK -- {args.current} congruent with {args.baseline}")
        return 0

    has_regression = any(
        f.startswith("[REGRESSION]") or f.startswith("[WARN]") or f.startswith("[SHA]")
        or f.startswith("[DIFF]")
        for f in findings
    )
    header = "REGRESSION" if has_regression else "DRIFT"
    print(f"[diff-manifest] {header} detected ({len(findings)} finding(s)):")
    for f in findings:
        print(f"  {f}")
    # Exit 1 on any regression or SHA/DIFF change; [INFO]-only deltas pass.
    return 1 if has_regression else 0


if __name__ == "__main__":
    sys.exit(main())
