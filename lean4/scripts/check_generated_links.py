#!/usr/bin/env python3
"""Phase I: broken-link checker for generated artifacts referenced by
the manuscript and supplemental material.

Scans every *.tex file in the author's manuscript source tree for
\\input{...}, \\include{...}, and Lean-source references like
\\texttt{KstarFormal/...}; confirms each referenced file exists on
disk.  Designed to catch drift when a Lean file gets renamed or a
generated artifact goes stale.  In a public clone the source tree
won't exist (author-only artifact) and the check skips cleanly.

Exit codes:
  0 -- all references resolve (or source tree absent = skip)
  1 -- one or more broken references
  2 -- invocation error
"""
from __future__ import annotations
import os
import re
import sys
from pathlib import Path

LEAN4_ROOT = Path(__file__).resolve().parent.parent
KSTAR_ROOT = LEAN4_ROOT.parent
SUBMISSION_ROOT = KSTAR_ROOT.parent
MANUSCRIPT_ROOT = Path(
    os.environ.get("KSTAR_MANUSCRIPT_ROOT", SUBMISSION_ROOT / "manuscript-source")
)

# `\input`/`\include` with optional bracketed args and a whitespace boundary
# so we don't match `\inputfoo{...}` (where `foo` is part of a longer command).
INPUT_RE = re.compile(r"\\(?:input|include)\b(?:\[[^\]]*\])?\{([^}]+)\}")
LEAN_REF_RE = re.compile(r"\\texttt\{(KstarFormal(?:[/\\][A-Za-z0-9_]+)+)(?:\.lean)?\}")


def find_tex_files() -> list[Path]:
    # Restrict to the companion manuscript tree; ignore backups and archives.
    skip = ("backup-", "ARCHIVED", "__pycache__")
    return [
        p
        for p in MANUSCRIPT_ROOT.rglob("*.tex")
        if not any(s in p.as_posix() for s in skip)
    ]


def resolve_input(rel: str, base: Path) -> Path | None:
    """Resolve a \\input{name} argument to a file.  LaTeX appends .tex
    if no extension is given."""
    candidates = [base / rel]
    if not rel.endswith(".tex"):
        candidates.append(base / (rel + ".tex"))
    for c in candidates:
        if c.is_file():
            return c
    return None


def resolve_lean_ref(ref: str) -> Path | None:
    # Normalise backslashes and append .lean
    normalized = ref.replace("\\", "/")
    if not normalized.endswith(".lean"):
        normalized += ".lean"
    p = LEAN4_ROOT / normalized
    return p if p.is_file() else None


def main() -> int:
    if not MANUSCRIPT_ROOT.is_dir():
        # Standalone kstar-verify clone (no manuscript tree) —
        # link-check is a no-op.  CI still runs it to catch drift when
        # the sibling tree IS present.
        print(f"[link-check] manuscript tree absent at {MANUSCRIPT_ROOT}; skipping")
        return 0
    tex_files = find_tex_files()
    if not tex_files:
        print(f"[link-check] no .tex files found under {MANUSCRIPT_ROOT}; skipping")
        return 0

    broken: list[tuple[Path, str, str]] = []
    ok_count = 0
    total = 0

    for tex in tex_files:
        try:
            text = tex.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        # \input / \include
        for m in INPUT_RE.finditer(text):
            total += 1
            arg = m.group(1)
            # Skip \input{<something>.bbl} or anything with whitespace
            # markers; LaTeX input typically single-line.
            resolved = resolve_input(arg, tex.parent)
            if resolved is None:
                broken.append((tex, arg, "\\input"))
            else:
                ok_count += 1
        # \texttt{KstarFormal/...} -- Lean source path references
        for m in LEAN_REF_RE.finditer(text):
            total += 1
            arg = m.group(1)
            resolved = resolve_lean_ref(arg)
            if resolved is None:
                broken.append((tex, arg, "Lean-ref"))
            else:
                ok_count += 1

    print(f"[link-check] scanned {len(tex_files)} .tex files, {total} references")
    if broken:
        print(f"[link-check] {len(broken)} BROKEN:")
        for tex, ref, kind in broken:
            try:
                rel_tex = tex.relative_to(SUBMISSION_ROOT)
            except ValueError:
                rel_tex = tex  # outside expected tree; print absolute
            print(f"  {kind}: {ref}  (in {rel_tex})")
        return 1
    print(f"[link-check] OK -- all {ok_count} references resolve")
    return 0


if __name__ == "__main__":
    sys.exit(main())
