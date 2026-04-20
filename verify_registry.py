#!/usr/bin/env python3
"""
Standalone registry verifier (no manuscript / proofs.tex dependency).

The registry (proofs_registry.yaml) is the SOLE artifact crossing the
trust boundary from the manuscript-side extraction pipeline into this
verification repo. Anchor substring
matching against the source TeX is enforced upstream by
verify_registry_anchors.py before the registry is signed and shipped
here. This script therefore checks ONLY what is verifiable from the
artifacts that ship with the repo:

  1. Registry schema  - every entry has the required fields, no unknown
                        fields, no duplicate IDs.
  2. Anchor sanity    - non-empty list, count in [2, 10], min length >= 8,
                        no duplicate anchors across the registry.
  3. depends_on graph - no dangling references, no cycles.
  4. Lean mapping     - lean4_theorem and lean4_tier1 names are present
                        in lean4/KstarFormal/**/*.lean as theorem/lemma
                        declarations.
  5. lean4_status     - is one of the recognised values.
  6. Provenance       - SHA-256 of the registry is printed.

No proofs.tex. No manuscript.tex. No --proofs-tex flag. The repo cannot
re-extract the proof text and must trust the upstream signed artifact.

Usage:
  python verify_registry.py
  python verify_registry.py --registry proofs_registry.yaml
"""
import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

try:
    import yaml
except ImportError:
    print("ERROR: pyyaml required (pip install pyyaml)")
    sys.exit(1)

ROOT = Path(__file__).parent.resolve()
LEAN_ROOT = ROOT / "lean4" / "KstarFormal"
KSTAR_LABELS_LEAN = LEAN_ROOT / "Combinatorics" / "GhzNonCoverage.lean"

REQUIRED = {"id", "label", "type", "lines", "claims", "claim_anchors", "depends_on"}
OPTIONAL = {
    "source",
    "known_caveats",
    "lean4_theorem",
    "lean4_tier1",
    "lean4_status",
    "lean4_axioms",
    # Phase D/G (2026-04-17): universal-n Lean wrappers.
    "lean4_theorem_universal",
    "lean4_universal_axioms",
    "lean4_universal_witnesses_n",
    "key_values",
    # Process / adversarial provenance fields (Item 7).
    # `mutation_killed`: bool — true iff at least one adversarial Lean
    #   mutation of the load-bearing theorem(s) has been generated and
    #   shown to fail to compile (i.e., the proof is not vacuously true).
    # `mutation_killed_count`: int — number of distinct mutations killed.
    "mutation_killed",
    "mutation_killed_count",
    # Extended provenance fields (added for fact-type entries and
    # enriched theorem entries).
    "key_value_roles",
    "proof_method",
    "verification_scope",
    "verified_by_tier",
    "cert_artifacts",
    "tier3_test",
    "note",
}
VALID_STATUS = {
    "proved",
    # Item 8: split `proved-from-axioms` into two sub-categories.
    # `proved-from-foundational` — proof relies only on Lean kernel /
    #     classical-logic axioms (propext, Quot.sound, Classical.choice).
    #     These are the universally accepted Mathlib foundations and do
    #     not represent paper-specific mathematical assumptions.
    # `proved-from-mathematical` — proof additionally invokes at least
    #     one mathematical axiom declared in KstarFormal.Axioms (e.g.,
    #     Uhlmann, Schatten, Parseval, Weyl). These are textbook math
    #     results assumed without re-proof and are the load-bearing
    #     mathematical trust base of the paper.
    "proved-from-foundational",
    "proved-from-mathematical",
    # Deprecated alias for backwards compatibility; equivalent to
    # `proved-from-mathematical` for verification purposes. New entries
    # should use one of the two finer-grained labels above.
    "proved-from-axioms",
    "informal",
    "pending",
    "tier1-pending",
    "not_required",
}
VALID_TYPES = {"lemma", "theorem", "proposition", "corollary", "fact", "equation"}

PASS = 0
FAIL = 0


def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  [PASS] {name}" + (f"  ({detail})" if detail else ""))
    else:
        FAIL += 1
        print(f"  [FAIL] {name}" + (f"  ({detail})" if detail else ""))


def _collect_lean_decls(lean_root: Path) -> set[str]:
    """Return the set of theorem/lemma names declared anywhere under
    lean_root, excluding the .lake/ build directory."""
    pat = re.compile(r"^\s*(?:theorem|lemma)\s+([A-Za-z_][\w']*)", re.MULTILINE)
    names: set[str] = set()
    if not lean_root.is_dir():
        return names
    for f in lean_root.rglob("*.lean"):
        if ".lake" in f.parts:
            continue
        try:
            text = f.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        names.update(pat.findall(text))
    return names


def _kstar_labels_canonical_sha(lean_file: Path) -> str | None:
    """Parse `kstar_labels_n4` out of the Lean source and compute the
    canonical SHA-256 over its JSON serialization. Returns None if the
    file or constant cannot be parsed."""
    if not lean_file.is_file():
        return None
    text = lean_file.read_text(encoding="utf-8", errors="replace")
    m = re.search(
        r"def\s+kstar_labels_n4\s*:\s*List\s*\(List Nat\)\s*:=\s*\[(.*?)\]\s*\n\s*\n",
        text,
        re.DOTALL,
    )
    if m is None:
        return None
    nums = [int(n) for n in re.findall(r"\d+", m.group(1))]
    if not nums or len(nums) % 4 != 0:
        return None
    canon = [nums[i : i + 4] for i in range(0, len(nums), 4)]
    return hashlib.sha256(
        json.dumps(canon, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _has_cycle(graph: dict[str, list[str]]) -> str | None:
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {n: WHITE for n in graph}

    def dfs(n: str) -> str | None:
        color[n] = GRAY
        for m in graph.get(n, []):
            if m not in color:
                continue
            if color[m] == GRAY:
                return f"{n} -> {m}"
            if color[m] == WHITE:
                r = dfs(m)
                if r:
                    return r
        color[n] = BLACK
        return None

    for n in graph:
        if color[n] == WHITE:
            r = dfs(n)
            if r:
                return r
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        default=str(ROOT / "proofs_registry.yaml"),
        help="Path to registry YAML (default: proofs_registry.yaml in repo root)",
    )
    args = parser.parse_args()

    registry_path = Path(args.registry).resolve()
    if not registry_path.is_file():
        print(f"ERROR: registry not found: {registry_path}")
        return 2

    reg_bytes = registry_path.read_bytes()
    reg_sha = hashlib.sha256(reg_bytes).hexdigest()
    registry = yaml.safe_load(reg_bytes.decode("utf-8"))

    statements = registry.get("statements") or []

    print("=" * 70)
    print("STANDALONE REGISTRY VERIFICATION (no proofs.tex required)")
    print("=" * 70)
    print(f"  registry:    {registry_path}")
    print(f"  reg SHA-256: {reg_sha}")
    print(f"  statements:  {len(statements)}")
    print(f"  lean root:   {LEAN_ROOT}")
    print()

    check("Registry has at least one statement", len(statements) > 0)

    # ----- 1. duplicate ids -----
    ids = [s["id"] for s in statements]
    dup_ids = [i for i, c in Counter(ids).items() if c > 1]
    check("No duplicate statement IDs", not dup_ids,
          f"duplicates: {dup_ids}" if dup_ids else "")

    id_set = set(ids)

    # ----- 2. per-statement schema -----
    print("\n--- schema ---")
    for s in statements:
        sid = s.get("id", "<missing>")
        keys = set(s)
        missing = REQUIRED - keys
        unknown = keys - REQUIRED - OPTIONAL
        check(f"{sid}: required fields present", not missing,
              f"missing {sorted(missing)}" if missing else "")
        check(f"{sid}: no unknown fields", not unknown,
              f"unknown {sorted(unknown)}" if unknown else "")
        if s.get("type") not in VALID_TYPES:
            check(f"{sid}: type is recognised", False,
                  f"got {s.get('type')!r}")
        lines_val = s.get("lines")
        if not isinstance(lines_val, list) or len(lines_val) != 2:
            check(f"{sid}: lines is [start, end]", False,
                  f"got {lines_val!r}")
        elif (not isinstance(lines_val[0], int) or not isinstance(lines_val[1], int)
              or lines_val[0] < 1 or lines_val[1] < lines_val[0]):
            check(f"{sid}: lines are positive ints with start <= end", False,
                  f"got {lines_val!r}")
        if not s.get("claims"):
            check(f"{sid}: claims is non-empty", False)
        if not s.get("label"):
            check(f"{sid}: label is non-empty", False)

    # ----- 3. anchor sanity -----
    print("\n--- anchors ---")
    all_anchors: list[tuple[str, str]] = []
    for s in statements:
        sid = s["id"]
        anchors = s.get("claim_anchors") or []
        check(f"{sid}: has anchors", len(anchors) > 0)
        # Equation-type entries legitimately have 1-2 anchors (the LaTeX
        # equation label + optional symbol); other statement types should
        # have 2-10 anchor phrases for drift detection.
        min_anchors = 1 if s.get("type") == "equation" else 2
        check(f"{sid}: anchor count in [{min_anchors}, 10]",
              min_anchors <= len(anchors) <= 10,
              f"got {len(anchors)}")
        for a in anchors:
            all_anchors.append((sid, a))
            if len(a) < 8:
                check(f"{sid}: anchor length >= 8", False, f"{a!r}")
    counts = Counter(a for _, a in all_anchors)
    # Anchors shared between a theorem/lemma and an `equation` type are
    # expected: the equation entry IS the shared expression, and the
    # theorem cites it.  Only flag truly spurious sharing.
    eq_ids = {s["id"] for s in statements if s.get("type") == "equation"}
    shared: dict[str, int] = {}
    for a, c in counts.items():
        if c <= 1:
            continue
        owners = {sid for sid, aa in all_anchors if aa == a}
        if owners & eq_ids:
            # At least one owner is an equation; legitimate cross-ref.
            continue
        shared[a] = c
    check("No spurious anchor sharing between non-equation statements",
          not shared, f"{len(shared)} shared" if shared else "")

    # ----- 4. depends_on graph -----
    print("\n--- depends_on ---")
    graph = {s["id"]: list(s.get("depends_on") or []) for s in statements}
    dangling: list[tuple[str, str]] = []
    for sid, deps in graph.items():
        for d in deps:
            if d not in id_set:
                dangling.append((sid, d))
    check("All depends_on targets exist in registry", not dangling,
          f"dangling: {dangling}" if dangling else "")
    cyc = _has_cycle(graph)
    check("depends_on graph is acyclic", cyc is None,
          f"cycle through {cyc}" if cyc else "")

    # ----- 5. Lean mapping -----
    print("\n--- lean4 ---")
    decls = _collect_lean_decls(LEAN_ROOT)
    check(f"Lean tree present ({len(decls)} theorem/lemma decls)", len(decls) > 0)
    for s in statements:
        sid = s["id"]
        for fld in ("lean4_theorem", "lean4_tier1"):
            name = s.get(fld)
            if not name:
                continue
            check(f"{sid}.{fld}={name} present in lean4/", name in decls)
        status = s.get("lean4_status")
        if status is not None:
            check(f"{sid}: lean4_status is recognised",
                  status in VALID_STATUS,
                  f"got {status!r}" if status not in VALID_STATUS else "")
            axioms = s.get("lean4_axioms")
            # Item 8: foundational/mathematical consistency check.
            if status == "proved-from-mathematical":
                check(
                    f"{sid}: proved-from-mathematical has non-empty lean4_axioms",
                    isinstance(axioms, list) and len(axioms) > 0,
                    f"got {axioms!r}",
                )
            elif status == "proved-from-foundational":
                check(
                    f"{sid}: proved-from-foundational has empty lean4_axioms",
                    isinstance(axioms, list) and len(axioms) == 0,
                    f"got {axioms!r}",
                )

    # ----- 6. mutation_killed coverage (Item 7) -----
    print("\n--- mutation_killed ---")
    mk_total = 0
    mk_true = 0
    mk_count_total = 0
    for s in statements:
        sid = s["id"]
        if "mutation_killed" not in s and "mutation_killed_count" not in s:
            continue
        mk_total += 1
        mk = s.get("mutation_killed")
        if mk is not None:
            check(f"{sid}: mutation_killed is bool",
                  isinstance(mk, bool),
                  f"got {type(mk).__name__}")
            if mk is True:
                mk_true += 1
        mkc = s.get("mutation_killed_count")
        if mkc is not None:
            check(f"{sid}: mutation_killed_count is non-negative int",
                  isinstance(mkc, int) and mkc >= 0,
                  f"got {mkc!r}")
            if isinstance(mkc, int) and mkc >= 0:
                mk_count_total += mkc
    if mk_total == 0:
        print("  (no entries declare mutation_killed yet — schema slot active)")
    else:
        print(f"  coverage: {mk_true}/{mk_total} entries with mutation_killed=true")
        print(f"  total mutations killed: {mk_count_total}")

    # ----- 7. K* Lean<->Python bridge SHA cross-check -----
    print("\n--- kstar bridge ---")
    bridge_entry = next(
        (s for s in statements if s.get("id") == "fact:kstar_python_lean_bridge"),
        None,
    )
    if bridge_entry is None:
        check("fact:kstar_python_lean_bridge entry present", False,
              "no bridge entry in registry")
    else:
        pinned = (bridge_entry.get("key_values") or {}).get("canonical_sha256")
        check("bridge entry has canonical_sha256 in key_values",
              isinstance(pinned, str) and len(pinned) == 64,
              f"got {pinned!r}")
        derived = _kstar_labels_canonical_sha(KSTAR_LABELS_LEAN)
        check(f"Lean kstar_labels_n4 parseable from {KSTAR_LABELS_LEAN.name}",
              derived is not None)
        if pinned and derived:
            check("registry-pinned SHA == SHA(Lean kstar_labels_n4)",
                  pinned == derived,
                  f"pinned={pinned[:16]}... derived={derived[:16]}...")

    # ----- summary -----
    print(f"\n{'=' * 70}")
    print(f"REGISTRY VERIFICATION: {PASS} passed, {FAIL} failed")
    print(f"Registry SHA-256: {reg_sha}")
    print("=" * 70)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
