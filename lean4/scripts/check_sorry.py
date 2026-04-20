#!/usr/bin/env python3
"""Sorry audit for Lean4 proof files.

Layer 1 files must be sorry-free (complete proofs).
Layer 2/3 files may have sorry (stubs for future work).

Usage:
  python check_sorry.py             # report all sorry
  python check_sorry.py --layer1    # exit 1 if Layer 1 has sorry
  python check_sorry.py --stats     # include theorem/lemma counts per file
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

LAYER1_FILES = [
    "KstarFormal/Defs.lean",
    "KstarFormal/Combinatorics/Krawtchouk.lean",
    "KstarFormal/Combinatorics/LatticeCount.lean",
    "KstarFormal/Combinatorics/GreedyRedist.lean",
    "KstarFormal/Combinatorics/WeightSaturation.lean",
    "KstarFormal/LinearAlgebra/GramMatrix.lean",
    "KstarFormal/LinearAlgebra/SpectralDecomp.lean",
    "KstarFormal/LinearAlgebra/Eigenvalues.lean",
    "KstarFormal/LinearAlgebra/Monotonicity.lean",
    "KstarFormal/Certified.lean",
]

LAYER23_FILES = [
    "KstarFormal/LinearAlgebra/PauliOrthogonality.lean",
    "KstarFormal/Quantum/PurityBound.lean",
    "KstarFormal/Quantum/FidelityDichotomy.lean",
    "KstarFormal/Quantum/BasinSeparation.lean",
    "KstarFormal/Probability/Hypergeometric.lean",
    "KstarFormal/Statements.lean",
]

# Match 'sorry' as a Lean keyword, not inside line comments (--)
SORRY_RE = re.compile(r"^(?!\s*--).*\bsorry\b", re.MULTILINE)

# Block comments: /- ... -/  (may span lines)
BLOCK_COMMENT_RE = re.compile(r"/-.*?-/", re.DOTALL)

# Theorem/lemma declarations
DECL_RE = re.compile(r"^(theorem|lemma|def)\s+\w+", re.MULTILINE)


def _strip_block_comments(text):
    """Remove /- ... -/ block comments (may be nested in practice,
    but single-level stripping catches 'sorry' in docstrings)."""
    return BLOCK_COMMENT_RE.sub("", text)


def count_sorry(filepath):
    """Count sorry occurrences outside of line comments and block comments."""
    if not filepath.is_file():
        return -1  # missing
    text = filepath.read_text(encoding="utf-8", errors="replace")
    text = _strip_block_comments(text)
    return len(SORRY_RE.findall(text))


def count_declarations(filepath):
    """Count theorem/lemma/def declarations in a file."""
    if not filepath.is_file():
        return 0
    text = filepath.read_text(encoding="utf-8", errors="replace")
    return len(DECL_RE.findall(text))


def main():
    layer1_only = "--layer1" in sys.argv
    show_stats = "--stats" in sys.argv

    print("=== Sorry Audit ===")
    print()

    # Layer 1
    print("--- Layer 1 (target: sorry-free) ---")
    l1_sorry = 0
    l1_results = []
    l1_decls = 0
    for rel in LAYER1_FILES:
        fp = ROOT / rel
        n = count_sorry(fp)
        decls = count_declarations(fp) if show_stats else 0
        l1_decls += decls
        suffix = f" ({decls} decls)" if show_stats else ""
        if n < 0:
            print(f"  [MISS]  {rel}: not found")
            l1_results.append(("MISS", rel))
        elif n > 0:
            print(f"  [SORRY] {rel}: {n}{suffix}")
            l1_sorry += n
            l1_results.append(("SORRY", rel))
        else:
            print(f"  [CLEAN] {rel}{suffix}")
            l1_results.append(("CLEAN", rel))

    print()

    # Layer 2/3
    print("--- Layer 2/3 (sorry expected) ---")
    l23_sorry = 0
    l23_decls = 0
    for rel in LAYER23_FILES:
        fp = ROOT / rel
        n = count_sorry(fp)
        decls = count_declarations(fp) if show_stats else 0
        l23_decls += decls
        suffix = f" ({decls} decls)" if show_stats else ""
        if n < 0:
            print(f"  [MISS]  {rel}: not found")
        elif n > 0:
            print(f"  [STUB]  {rel}: {n}{suffix}")
            l23_sorry += n
        else:
            print(f"  [CLEAN] {rel}{suffix}")

    print()
    print("--- Summary ---")
    print(f"  Layer 1 sorry: {l1_sorry}")
    print(f"  Layer 2/3 sorry: {l23_sorry}")
    print(f"  Total sorry: {l1_sorry + l23_sorry}")
    if show_stats:
        print(f"  Layer 1 declarations: {l1_decls}")
        print(f"  Layer 2/3 declarations: {l23_decls}")
        print(f"  Total declarations: {l1_decls + l23_decls}")

    if layer1_only:
        if l1_sorry > 0:
            print(f"\nFAIL: Layer 1 has {l1_sorry} sorry")
            sys.exit(1)
        else:
            msg = f"PASS: Layer 1 is sorry-free"
            if show_stats:
                msg += f" ({l1_decls} theorems/lemmas verified)"
            print(f"\n{msg}")
            sys.exit(0)


if __name__ == "__main__":
    main()
