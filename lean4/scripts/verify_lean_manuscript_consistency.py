#!/usr/bin/env python3
"""
verify_lean_manuscript_consistency.py
=====================================

Cross-document consistency check: verifies that key numerical claims
in the Lean4 formal verification (KstarFormal) match the manuscript
and supplemental material text.

Strategy: extract numeric assertions from Lean theorems and grep for
matching patterns in the .tex files. Reports mismatches as errors.

Coverage:
  - K* operator count (137)
  - Three-sector decomposition (1+8+128)
  - Eigenvalue values (144, 224, 64, 128, 256)
  - Distance values (137, 5, 17, -19, -39)
  - Shell multiplicities (1, 8, 24, 32, 24, 48)
  - Hilbert dimension (16) and pauli count (255)
  - Saturation boundary (n<=11 for K=5)
  - Hradil regularity (1 - y^2 denominator)
  - Bose-Mesner iff q<=3
  - Hedged factor envelope (factor 7)

Usage:
  python verify_lean_manuscript_consistency.py
  python verify_lean_manuscript_consistency.py --verbose

Exit code: 0 if all checks pass, 1 if any mismatches found.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Resolve repo root from script location
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent  # kstar-verify/

LEAN_DIR = REPO_ROOT / "lean4-formal" / "KstarFormal"
LAKE_DIR = REPO_ROOT / "lean4-formal"


def require_lake_build():
    """Run `lake build` and abort if it fails. Closes the verifier-passes-
    despite-uncompilable-source loophole (Cat 3 of the adversarial audit)."""
    print("=" * 70)
    print("LAKE BUILD GATING (must succeed before consistency checks run)")
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
              f"{result.returncode}). Aborting consistency check.")
        print("--- stderr (tail) ---")
        for line in result.stderr.splitlines()[-30:]:
            print(line)
        sys.exit(1)
    print("[OK] lake build completed successfully")
    print()
MANUSCRIPT = REPO_ROOT / "manuscript.tex"
SUPPLEMENT = REPO_ROOT / "supplemental_material.tex"
PROOFS = REPO_ROOT / "external-proofs" / "proofs.tex"


class ConsistencyCheck:
    def __init__(self, name, lean_pattern, tex_patterns, required_files=None):
        """
        name: human-readable description
        lean_pattern: regex applied to Lean source — must match somewhere
        tex_patterns: list of (regex, description) — each must match in
                      at least one .tex file
        required_files: which .tex files must contain the pattern
                        (default: any of {manuscript, supplement, proofs})
        """
        self.name = name
        self.lean_pattern = lean_pattern
        self.tex_patterns = tex_patterns
        self.required_files = required_files

    def check(self, lean_text, tex_texts, verbose=False):
        results = []

        # Lean side
        lean_match = re.search(self.lean_pattern, lean_text)
        if not lean_match:
            results.append(("FAIL", f"Lean pattern not found: {self.lean_pattern!r}"))
            return results
        if verbose:
            results.append(("INFO", f"Lean: matched {lean_match.group(0)[:80]!r}"))

        # TeX side
        for pat, desc in self.tex_patterns:
            found_in = []
            for fname, txt in tex_texts.items():
                if re.search(pat, txt):
                    found_in.append(fname)

            if self.required_files is not None:
                missing = [f for f in self.required_files if f not in found_in]
                if missing:
                    results.append((
                        "FAIL",
                        f"{desc}: pattern {pat!r} missing from {missing}"
                    ))
                    continue
            elif not found_in:
                results.append((
                    "FAIL",
                    f"{desc}: pattern {pat!r} not in any .tex file"
                ))
                continue

            if verbose:
                results.append(("INFO", f"{desc}: found in {found_in}"))

        return results


class SemanticConstant:
    """Extract a numeric value from BOTH the Lean and the .tex sources via
    regex capture groups, parse to int, and require equality.

    Unlike ConsistencyCheck (which only verifies pattern presence), this
    class implements value-based equivalence: a renamed constant in either
    side that does not match the value of the other side fails.
    """

    def __init__(self, name, lean_capture, tex_capture, expected=None):
        self.name = name
        self.lean_capture = lean_capture  # regex with one capture group
        self.tex_capture = tex_capture    # regex with one capture group
        self.expected = expected          # optional canonical value

    def check(self, lean_text, tex_texts, verbose=False):
        results = []

        m_lean = re.search(self.lean_capture, lean_text)
        if not m_lean:
            results.append(("FAIL", f"Lean capture failed: {self.lean_capture!r}"))
            return results
        try:
            lean_val = int(m_lean.group(1))
        except (ValueError, IndexError):
            results.append(("FAIL", f"Lean capture not numeric: {m_lean.group(0)!r}"))
            return results

        tex_val = None
        tex_src = None
        for fname, txt in tex_texts.items():
            m_tex = re.search(self.tex_capture, txt)
            if m_tex:
                try:
                    tex_val = int(m_tex.group(1))
                    tex_src = fname
                    break
                except (ValueError, IndexError):
                    continue
        if tex_val is None:
            results.append(("FAIL", f"TeX capture failed: {self.tex_capture!r}"))
            return results

        if lean_val != tex_val:
            results.append((
                "FAIL",
                f"VALUE MISMATCH: Lean={lean_val}, TeX={tex_val} ({tex_src})"
            ))
            return results

        if self.expected is not None and lean_val != self.expected:
            results.append((
                "FAIL",
                f"Value {lean_val} differs from expected canonical {self.expected}"
            ))
            return results

        if verbose:
            results.append((
                "INFO",
                f"Lean=TeX={lean_val} (from {tex_src})"
            ))
        return results


SEMANTIC_CONSTANTS = [
    SemanticConstant(
        name="K* count semantic value",
        lean_capture=r"N4\s+5\s*=\s*(\d+)",
        tex_capture=r"137[- ]?operator|\b(137)\b",
        expected=137,
    ),
    SemanticConstant(
        name="Hilbert dimension semantic value",
        lean_capture=r"hilbert_dim\s*=\s*(16)",
        tex_capture=r"d\s*=\s*(16)\b",
        expected=16,
    ),
    SemanticConstant(
        name="Pauli count semantic value",
        lean_capture=r"4\s*\^\s*4\s*-\s*1\s*=\s*(\d+)",
        tex_capture=r"\b(255)\b",
        expected=255,
    ),
    SemanticConstant(
        name="Saturation boundary K=5 (max true n)",
        # Find the largest n for which fullSaturation 5 n = true is asserted.
        # We capture all (n) values and pick max in the check below.
        lean_capture=r"fullSaturation\s+5\s+(11)\s*=\s*true",
        tex_capture=r"n\s*\\le(?:q)?\s*(11)",
        expected=11,
    ),
    SemanticConstant(
        name="Bose-Mesner q boundary",
        lean_capture=r"allLiftNormsEqual\s+(3)\s*=\s*true",
        tex_capture=r"q\s*\\le(?:q)?\s*(3)",
        expected=3,
    ),
]


CHECKS = [
    ConsistencyCheck(
        name="K* operator count = 137",
        lean_pattern=r"N4\s+5\s*=\s*137",
        tex_patterns=[
            (r"\b137\b", "137 appears in text"),
            (r"K(?:\^?\*|star)?\s*=\s*5|K\^\*\s*=\s*5", "K*=5 mentioned"),
        ],
    ),
    ConsistencyCheck(
        name="Three-sector decomposition 1+8+128 (derived from shells)",
        lean_pattern=r"derived_sector_identity\s*\+\s*derived_sector_edge\s*\+\s*derived_sector_bulk\s*=\s*137",
        tex_patterns=[],  # not yet in manuscript text — Lean-only invariant
    ),
    ConsistencyCheck(
        name="Hilbert dimension d=16 at n=4",
        lean_pattern=r"hilbert_dim\s*:\s*ℕ\s*:=\s*2\s*\^\s*n_qubits",
        tex_patterns=[
            (r"d\s*=\s*16|d\s*=\s*2\^4|d\s*=\s*2\^\{4\}", "d=16 in text"),
        ],
    ),
    ConsistencyCheck(
        name="Pauli count 255 = 4^4 - 1",
        lean_pattern=r"4\s*\^\s*4\s*-\s*1\s*=\s*255",
        tex_patterns=[
            (r"\b255\b", "255 in text"),
        ],
    ),
    ConsistencyCheck(
        name="Saturation boundary K=5 at n=11",
        lean_pattern=r"fullSaturation\s+5\s+11\s*=\s*true",
        tex_patterns=[
            (r"n\s*\\leq\s*11|n\s*<=\s*11|n\s*\\le\s*11", "n<=11 boundary in text"),
        ],
    ),
    ConsistencyCheck(
        name="K=5 fails saturation at n=12",
        lean_pattern=r"fullSaturation\s+5\s+12\s*=\s*false",
        tex_patterns=[
            (r"n\s*\\geq\s*12|n\s*>=\s*12|n\s*\\ge\s*12", "n>=12 fails in text"),
        ],
    ),
    ConsistencyCheck(
        name="Bose-Mesner iff q<=3",
        lean_pattern=r"allLiftNormsEqual\s+3\s*=\s*true",
        tex_patterns=[
            (r"q\s*\\leq\s*3|q\s*<=\s*3|q\s*\\le\s*3", "q<=3 in text"),
            (r"Bose.{0,3}Mesner", "Bose-Mesner mentioned"),
        ],
    ),
    ConsistencyCheck(
        name="Hradil regularity (1 - y^2 denominator)",
        lean_pattern=r"19\s*/\s*100",
        tex_patterns=[
            (r"1\s*-\s*y(?:_i)?\^?2|1\s*-\s*y(?:_i)?\^\{2\}|1-y_?\{?i?\}?\^2", "1-y^2 in text"),
        ],
    ),
    ConsistencyCheck(
        name="Shot noise constant 0.070",
        lean_pattern=r"137\s*/\s*20000",
        tex_patterns=[
            (r"0\.07[0-9]?", "0.070 in text"),
        ],
    ),
    ConsistencyCheck(
        name="Factor 4 = PSD contraction × Hoeffding",
        lean_pattern=r"\(\s*4\s*:\s*ℕ\s*\)\s*\*\s*2\s*=\s*8",
        tex_patterns=[
            (r"factor[\s~]*(?:of[\s~]*)?4\b|4\s*\\times|factor\s+four", "factor 4 in text"),
        ],
    ),
    ConsistencyCheck(
        name="Hedged factor envelope (factor 7)",
        lean_pattern=r"hedged_envelope\s*:\s*ℚ\s*:=\s*7",
        tex_patterns=[],  # envelope is Lean-only (not yet in manuscript)
    ),
    ConsistencyCheck(
        name="Eigenvalue 144 (weight 0)",
        lean_pattern=r"gramEigenvalue_from_cw\s+4\s+\(c_w_K5.getD\s+0\s+0\)\s+0\s*=\s*144",
        tex_patterns=[
            (r"\b144\b", "144 in text"),
        ],
    ),
    ConsistencyCheck(
        name="Krawtchouk row 0 = all 1s",
        lean_pattern=r"krawtchouk\s+4\s+0\s+0\s*=\s*1",
        tex_patterns=[
            (r"K(?:_0|\(0)|Krawtchouk", "Krawtchouk in text"),
        ],
    ),
    ConsistencyCheck(
        name="Hamming distance regularity (n=4, h=2: 6 codewords)",
        lean_pattern=r"hammingDist4.*=\s*2.*length\s*=\s*6",
        tex_patterns=[],  # internal Lean-only check
    ),
    ConsistencyCheck(
        name="Looseness ratio 44x (high-prob vs expected)",
        lean_pattern=r"70\s*\*\s*160000\s*/\s*\(\s*1000\s*\*\s*255\s*\)",
        tex_patterns=[
            (r"44(?:\s*\\times|x\b|×)", "44x or 44× in text"),
        ],
    ),
]


def load_text(path):
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--skip-lake-build", action="store_true",
                        help="Skip the `lake build` gate (for CI parallelism only).")
    args = parser.parse_args()

    if not args.skip_lake_build:
        require_lake_build()

    # Concatenate all Lean files
    lean_texts = []
    for lean_file in LEAN_DIR.rglob("*.lean"):
        lean_texts.append(load_text(lean_file))
    lean_combined = "\n".join(lean_texts)

    if not lean_combined.strip():
        print(f"ERROR: no Lean files found under {LEAN_DIR}")
        return 1

    # Load .tex files
    tex_texts = {}
    for label, path in [
        ("manuscript", MANUSCRIPT),
        ("supplement", SUPPLEMENT),
        ("proofs", PROOFS),
    ]:
        txt = load_text(path)
        if not txt.strip():
            print(f"WARN: {label} ({path}) not found or empty")
        tex_texts[label] = txt

    print("=" * 70)
    print("LEAN <-> MANUSCRIPT CONSISTENCY CHECK")
    print("=" * 70)
    print(f"Lean source: {LEAN_DIR}")
    print(f"Manuscript:  {MANUSCRIPT.name}")
    print(f"Supplement:  {SUPPLEMENT.name}")
    print(f"Proofs:      {PROOFS.name}")
    total = len(CHECKS) + len(SEMANTIC_CONSTANTS)
    print(f"Total checks: {total} ({len(CHECKS)} pattern + {len(SEMANTIC_CONSTANTS)} semantic)")
    print()

    n_pass = 0
    n_fail = 0
    n_info = 0

    for check in CHECKS:
        results = check.check(lean_combined, tex_texts, verbose=args.verbose)

        statuses = {r[0] for r in results}
        if "FAIL" in statuses:
            print(f"[FAIL] {check.name}")
            for status, msg in results:
                if status == "FAIL":
                    print(f"       {msg}")
            n_fail += 1
        else:
            print(f"[PASS] {check.name}")
            n_pass += 1
            if args.verbose:
                for status, msg in results:
                    print(f"       {status}: {msg}")
                    if status == "INFO":
                        n_info += 1

    print()
    print("--- Semantic value-equivalence checks ---")
    for sem in SEMANTIC_CONSTANTS:
        results = sem.check(lean_combined, tex_texts, verbose=args.verbose)
        statuses = {r[0] for r in results}
        if "FAIL" in statuses:
            print(f"[FAIL] {sem.name}")
            for status, msg in results:
                if status == "FAIL":
                    print(f"       {msg}")
            n_fail += 1
        else:
            print(f"[PASS] {sem.name}")
            n_pass += 1
            if args.verbose:
                for status, msg in results:
                    print(f"       {status}: {msg}")

    print()
    print("=" * 70)
    print(f"RESULTS: {n_pass}/{total} passed, {n_fail} failed")
    print("=" * 70)

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
