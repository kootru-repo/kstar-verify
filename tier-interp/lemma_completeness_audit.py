r"""
lemma_completeness_audit.py - tex <-> Lean label coverage audit.

Walks the manuscript .tex files, extracts every \label{lem|thm|prop|cor:...},
and confirms each one is either (a) mapped to a concrete Lean theorem name
that actually exists in the lean4/ tree, or (b) explicitly classified as
informal (not formalizable in Tier 1 scope, with a reason).

Reports three failure modes:
  1. UNCLASSIFIED: tex label has no manifest entry.
  2. ROTTED: manifest claims label X maps to Lean theorem Y, but Y is no
     longer present in the Lean source.
  3. STALE: manifest references a tex label that no longer exists in the
     current .tex sources.

Path resolution:
  - Defaults BASE to `<repo>/manuscript-tex` and LEAN_ROOT to `<repo>/lean4`,
    relative to this script. Override with KSTAR_TEX_BASE and KSTAR_LEAN_ROOT.
  - If the default tex base does not exist, the script reports STALE/missing
    tex inputs as a soft warning rather than a hard failure (the public repo
    may ship without manuscript sources).

Expected output: 0 unclassified, 0 rotted, 0 stale -> exit 0.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TEX_BASE = (SCRIPT_DIR.parent / "manuscript-tex").resolve()
DEFAULT_LEAN_ROOT = (SCRIPT_DIR.parent / "lean4").resolve()

BASE = Path(os.environ.get("KSTAR_TEX_BASE", str(DEFAULT_TEX_BASE)))
LEAN_ROOT = Path(os.environ.get("KSTAR_LEAN_ROOT", str(DEFAULT_LEAN_ROOT)))

TEX_FILES = [
    BASE / "manuscript.tex",
    BASE / "supplemental_material.tex",
    BASE / "external-proofs" / "proofs.tex",
]

LABEL = re.compile(r"\\label\{((?:lem|thm|prop|cor|app):[^}]+)\}")

# Manifest: tex-label -> (Lean theorem name | None, status, note)
#   status in {"formalized", "informal", "tier1-pending"}
MANIFEST: dict[str, tuple[str | None, str, str]] = {
    "lem:hessian": (
        "lem1_hessian_diagonal",
        "formalized",
        "Statements.lean (Layer 1): Hessian diagonality + h_P formula",
    ),
    "thm:basin": (
        "thm1_iii_full_finite_sample_bound",
        "formalized",
        "BasinSeparation.lean",
    ),
    "cor:approx_local": (
        "thm1_iii_concrete_n4_pure_local",
        "formalized",
        "BasinSeparation.lean concrete instance",
    ),
    "prop:purity_main": (
        "lem2_purity_bound",
        "formalized",
        "Statements.lean (Layer 1): purity bound on positivity excursion",
    ),
    "lem:monotone": (
        "lem3_eigenvalue_monotonicity_n4",
        "formalized",
        "Statements.lean (Layer 1): eigenvalue monotonicity at n=4",
    ),
    "thm:spectral-char": (
        "bose_mesner_combinatorial_iff",
        "formalized",
        "Tier 1A: q-ary iff characterization (QaryGram.lean)",
    ),
    "prop:coupon": (
        "lem4_hypergeometric_bound",
        "formalized",
        "Statements.lean (Layer 1): hypergeometric / coupon-collector bound",
    ),
    "cor:lower_bound": (
        "cor2_operator_lower_bound_n4",
        "formalized",
        "Statements.lean (Layer 1): operator lower bound at n=4",
    ),
    "thm:asymptotic": (
        "thm3_asymptotic_separation",
        "formalized",
        "Statements.lean (Layer 1): asymptotic K* separation",
    ),
    "prop:spectral_q_main": (
        "bose_mesner_combinatorial_iff",
        "formalized",
        "Tier 1A: same iff characterization, manuscript-side label",
    ),
    "prop:spectral_q": (
        "bose_mesner_combinatorial_iff",
        "formalized",
        "Tier 1A iff characterization (SM-side label)",
    ),
    "cor:support": (
        "bose_mesner_combinatorial_iff",
        "formalized",
        "Tier 1A: support-completeness follows from iff characterization (QaryGram.lean)",
    ),
    "prop:spec-complete": (
        None,
        "informal",
        "Spectral completeness; covered by spectral_char_iff_q_le_3 in spirit",
    ),
    "prop:mass-threshold": (
        "weight12_saturation_K5_all_n",
        "formalized",
        "Tier 1B: WeightSatAllN.lean",
    ),
    "app:completeness": (
        "lem6_support_completeness_n4",
        "formalized",
        "Lemma 6 (support-completeness criterion); Statements.lean",
    ),
    "app:pauli-set": (
        None,
        "informal",
        "Appendix section label for the per-weight operator budget table; not a formal statement",
    ),
    "prop:psd_contraction": (
        None,
        "informal",
        "PSD contraction property; standard linear-algebra",
    ),
}


def extract_tex_labels() -> dict[str, list[Path]]:
    out: dict[str, list[Path]] = {}
    for f in TEX_FILES:
        if not f.exists():
            print(f"  WARNING: {f} not found", file=sys.stderr)
            continue
        text = f.read_text(encoding="utf-8")
        for m in LABEL.finditer(text):
            out.setdefault(m.group(1), []).append(f)
    return out


def find_lean_theorem(name: str) -> Path | None:
    pat = re.compile(rf"^\s*(?:theorem|lemma)\s+{re.escape(name)}\b", re.MULTILINE)
    for f in LEAN_ROOT.rglob("*.lean"):
        try:
            if pat.search(f.read_text(encoding="utf-8")):
                return f
        except Exception:
            continue
    return None


def main() -> int:
    print("=" * 72)
    print("LEMMA COMPLETENESS AUDIT - tex labels vs Lean theorems")
    print("=" * 72)
    print(f"  TEX_BASE:  {BASE}")
    print(f"  LEAN_ROOT: {LEAN_ROOT}")

    if not LEAN_ROOT.exists():
        print(f"  ERROR: LEAN_ROOT does not exist. Set KSTAR_LEAN_ROOT or run "
              f"from inside the repo.")
        return 2

    tex_present = BASE.exists() and any(f.exists() for f in TEX_FILES)
    if not tex_present:
        print(f"  NOTE: tex sources not present at {BASE}. Verifying Lean side "
              f"only (manifest -> Lean theorems).")

    tex_labels = extract_tex_labels() if tex_present else {}
    if tex_present:
        print(f"  Found {len(tex_labels)} unique labels across "
              f"{len(TEX_FILES)} .tex files")

    unclassified: list[str] = []
    rotted: list[tuple[str, str]] = []
    stale: list[str] = []
    formalized: list[tuple[str, str, Path]] = []
    informal: list[tuple[str, str]] = []
    pending: list[tuple[str, str]] = []

    # Verify manifest entries against Lean source.
    for label, (lean_name, status, note) in MANIFEST.items():
        if status == "formalized":
            assert lean_name is not None
            path = find_lean_theorem(lean_name)
            if path is None:
                rotted.append((label, lean_name))
            else:
                formalized.append((label, lean_name, path))
        elif status == "informal":
            informal.append((label, note))
        elif status == "tier1-pending":
            assert lean_name is not None
            path = find_lean_theorem(lean_name)
            if path is not None:
                formalized.append((label, lean_name, path))
            else:
                pending.append((label, note))

    # Cross-check tex side if available.
    if tex_present:
        for label in sorted(tex_labels):
            if label not in MANIFEST:
                unclassified.append(label)
        for label in MANIFEST:
            if label not in tex_labels:
                stale.append(label)

    print()
    print(f"  formalized:    {len(formalized)}")
    print(f"  informal:      {len(informal)}")
    print(f"  tier1-pending: {len(pending)}")
    print(f"  UNCLASSIFIED:  {len(unclassified)}")
    print(f"  ROTTED:        {len(rotted)}")
    print(f"  STALE:         {len(stale)}")

    if formalized:
        print()
        print("FORMALIZED (verified Lean theorem exists):")
        for label, name, path in formalized:
            try:
                rel = path.relative_to(LEAN_ROOT)
            except ValueError:
                rel = path
            print(f"  {label}  ->  {name}   [{rel}]")

    if unclassified:
        print()
        print("UNCLASSIFIED (no manifest entry - must classify):")
        for label in unclassified:
            print(f"  {label}")

    if rotted:
        print()
        print("ROTTED (manifest -> Lean theorem missing):")
        for label, name in rotted:
            print(f"  {label}  ->  {name}  (NOT FOUND)")

    if stale:
        print()
        print("STALE (manifest entry whose tex label no longer exists):")
        for label in stale:
            print(f"  {label}")

    print("=" * 72)
    # If tex sources are not present, only ROTTED is a hard failure.
    if not tex_present:
        return 0 if not rotted else 1
    return 0 if not (unclassified or rotted or stale) else 1


if __name__ == "__main__":
    sys.exit(main())
