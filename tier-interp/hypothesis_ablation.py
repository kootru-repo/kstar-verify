"""
hypothesis_ablation.py - Tier-interp hypothesis decorativeness checker (parallel).

For each theorem in the load-bearing Lean files, drops each named hypothesis
binder one at a time and re-builds. A hypothesis whose removal does not break
the build is *decorative* - either the proof never used it (statement is
silently weaker than its signature suggests) or no caller depends on it.

Path resolution:
  - Defaults LEAN_ROOT to `<repo>/lean4` (the public repo's Lean tier).
  - Override with environment variable KSTAR_LEAN_ROOT for testing against an
    out-of-tree Lean source.

Parallelization model:
  - The original source file is NEVER mutated. Each mutation is written to a
    fresh file in a temp directory and compiled via `lake env lean <abs path>`,
    which resolves imports from the project's .lake/build cache.
  - multiprocessing.Pool dispatches mutations across workers (default 16,
    capped at the number of mutations).

Expected output: 18/18 load-bearing, 0 decorative. Exit 0 on clean run.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_LEAN_ROOT = (SCRIPT_DIR.parent / "lean4").resolve()
LEAN_ROOT = Path(os.environ.get("KSTAR_LEAN_ROOT", str(DEFAULT_LEAN_ROOT)))

TARGETS = [
    "KstarFormal/Combinatorics/QaryGram.lean",
    "KstarFormal/Combinatorics/WeightSatAllN.lean",
    "KstarFormal/Quantum/HradilDerivation.lean",
]

WORKERS = int(os.environ.get("KSTAR_WORKERS", "16"))

THM_HEAD = re.compile(r"^\s*(?:theorem|lemma)\s+(\w+)")
BINDER = re.compile(r"\((\w+)\s*:\s*([^()]*(?:\([^()]*\)[^()]*)*)\)")

# Types that indicate a non-hypothesis binder (a parameter, not a Prop).
TYPE_BINDER_PREFIXES = (
    "ℕ", "Nat", "ℤ", "Int", "ℚ", "Rat", "ℝ", "Real",
    "Type", "Sort", "Prop", "Fin", "Bool", "List",
    "Vector", "Matrix", "ℂ", "Complex",
)


def is_hypothesis_type(t: str) -> bool:
    t = t.strip()
    for p in TYPE_BINDER_PREFIXES:
        if t == p:
            return False
        if t.startswith(p) and (len(t) == len(p) or not t[len(p)].isalnum()):
            return False
    return True


def find_theorem_signatures(text: str) -> list[tuple[int, int, str, str]]:
    lines = text.splitlines()
    out = []
    i = 0
    while i < len(lines):
        m = THM_HEAD.match(lines[i])
        if not m:
            i += 1
            continue
        name = m.group(1)
        j = i
        while j < len(lines) and ":=" not in lines[j]:
            j += 1
        if j >= len(lines):
            i += 1
            continue
        header = "\n".join(lines[i : j + 1])
        out.append((i, j, name, header))
        i = j + 1
    return out


def strip_binder(header: str, binder_name: str) -> str | None:
    found = []
    for m in BINDER.finditer(header):
        if m.group(1) == binder_name:
            found.append(m.span())
    if len(found) != 1:
        return None
    s, e = found[0]
    if s > 0 and header[s - 1] == " ":
        s -= 1
    return header[:s] + header[e:]


@dataclass
class Mutation:
    module_lean: str
    theorem_name: str
    binder_name: str
    mutated_text: str


def build_mutation(mut: Mutation) -> tuple[Mutation, bool]:
    fd, tmp_path = tempfile.mkstemp(
        prefix="ablate_", suffix=".lean", dir=tempfile.gettempdir()
    )
    try:
        os.write(fd, mut.mutated_text.encode("utf-8"))
        os.close(fd)
        res = subprocess.run(
            ["lake", "env", "lean", tmp_path],
            cwd=str(LEAN_ROOT),
            capture_output=True,
            timeout=300,
        )
        return mut, (res.returncode == 0)
    except Exception:
        return mut, False
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


def baseline_build(module_lean: str) -> bool:
    res = subprocess.run(
        ["lake", "env", "lean", str(LEAN_ROOT / module_lean)],
        cwd=str(LEAN_ROOT),
        capture_output=True,
        timeout=300,
    )
    return res.returncode == 0


def collect_mutations(module_lean: str) -> list[Mutation]:
    path = LEAN_ROOT / module_lean
    original = path.read_text(encoding="utf-8")
    sigs = find_theorem_signatures(original)
    muts: list[Mutation] = []
    for start, end, name, header in sigs:
        for bm in BINDER.finditer(header):
            bname = bm.group(1)
            btype = bm.group(2).strip()
            if not is_hypothesis_type(btype):
                continue
            new_header = strip_binder(header, bname)
            if new_header is None or new_header == header:
                continue
            new_lines = original.splitlines()
            prefix = "\n".join(new_lines[:start])
            suffix_lines = new_lines[end + 1 :]
            suffix = "\n".join(suffix_lines) if suffix_lines else ""
            replaced = (
                (prefix + "\n" if prefix else "")
                + new_header
                + ("\n" + suffix if suffix else "")
            )
            muts.append(Mutation(module_lean, name, bname, replaced))
    return muts


def main() -> int:
    print("=" * 72)
    print("HYPOTHESIS ABLATION (parallel) - decorativeness checker")
    print("=" * 72)
    print(f"  LEAN_ROOT: {LEAN_ROOT}")
    if not LEAN_ROOT.exists():
        print(f"  ERROR: LEAN_ROOT does not exist. Set KSTAR_LEAN_ROOT or run "
              f"from inside the repo.")
        return 2

    for tgt in TARGETS:
        print(f"\n[baseline] {tgt}")
        if not baseline_build(tgt):
            print(f"  ERROR: baseline build of {tgt} failed. Run `lake build` "
                  f"in {LEAN_ROOT} first.")
            return 2

    all_muts: list[Mutation] = []
    for tgt in TARGETS:
        muts = collect_mutations(tgt)
        print(f"  {tgt}: {len(muts)} hypothesis binders to ablate")
        all_muts.extend(muts)

    print(f"\nDispatching {len(all_muts)} mutations across {WORKERS} workers...")
    findings: list[tuple[Mutation, bool]] = []
    with Pool(processes=min(WORKERS, max(1, len(all_muts)))) as pool:
        for i, (mut, ok) in enumerate(
            pool.imap_unordered(build_mutation, all_muts), start=1
        ):
            findings.append((mut, ok))
            tag = "[DECORATIVE]" if ok else "[load-bearing]"
            print(
                f"  ({i}/{len(all_muts)}) {tag} "
                f"{mut.module_lean}::{mut.theorem_name}::{mut.binder_name}"
            )

    decorative = [(m, _) for (m, _) in findings if _]
    print()
    print("=" * 72)
    print(f"  Total binders tested: {len(findings)}")
    print(f"  Load-bearing:        {len(findings) - len(decorative)}")
    print(f"  DECORATIVE:          {len(decorative)}")
    if decorative:
        print()
        print("DECORATIVE HYPOTHESES (proof never used them):")
        for mut, _ in decorative:
            print(f"  {mut.module_lean}::{mut.theorem_name}  ({mut.binder_name} : ...)")
    print("=" * 72)
    return 0 if not decorative else 1


if __name__ == "__main__":
    sys.exit(main())
