"""
constant_fuzzing.py - Tier-interp numeric-literal fuzzer (parallel).

For each integer literal in the *signature* of a Tier 1 theorem (between
`theorem`/`lemma` and `:=`), perturb it by +/-1 and re-build. The expected
outcome is that the build breaks: a surviving perturbation means the theorem
statement is independent of that literal, which usually indicates a typo, a
copy-paste defect, or an over-general statement.

Path resolution:
  - Defaults LEAN_ROOT to `<repo>/lean4` (the public repo's Lean tier).
  - Override with environment variable KSTAR_LEAN_ROOT.

Scope:
  - Only literals <= 99 are mutated. Larger literals (e.g., 4096-element
    masks, 65535-row matrix sizes) are skipped because +/-1 perturbations
    of those typically break unrelated invariants and produce noise rather
    than meaningful independence findings. The Tier 1 statements use only
    small literals (n, q, weights, K*) that fall within this range.

Allowlist:
  - `free_allowlist.json` (next to this script) records FREE findings that
    are *intrinsic* to the lift-norm function symmetries (reflection
    a <-> q-a and saturation q >= 2a => lns = a^2). These are not defects;
    they are mathematical facts about the underlying combinatorial object.
  - On a clean run, the set of FREE findings must equal the allowlist
    exactly. Drift in either direction is reported and exits non-zero.
  - Expected output: 69 FREE findings, all in allowlist => exit 0.

Parallelization model:
  - The original source file is NEVER mutated. Each mutation is written to a
    fresh file in a temp directory and compiled via `lake env lean <abs path>`,
    which resolves imports from the project's .lake/build cache.
  - multiprocessing.Pool dispatches mutations across workers (default 16).
"""
from __future__ import annotations

import json
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
ALLOWLIST_PATH = SCRIPT_DIR / "free_allowlist.json"

TARGETS = [
    "KstarFormal/Combinatorics/QaryGram.lean",
    "KstarFormal/Combinatorics/WeightSatAllN.lean",
    "KstarFormal/Quantum/HradilDerivation.lean",
]

WORKERS = int(os.environ.get("KSTAR_WORKERS", "16"))

THM_HEAD = re.compile(r"^[ \t]*(?:theorem|lemma)\s+(\w+)")
LITERAL = re.compile(r"(?<![\w.])(\d+)(?![\w.])")
SKIP_CONTEXT = re.compile(
    r"(Fin|Sort|Type|Sort\*|max|getD|get!|get\?|nthLe|take|drop)\s*$"
)
SKIP_AFTER_INDEX = re.compile(r"(getD)\s+\d+\s*$")


def find_signatures(text: str) -> list[tuple[int, int, str]]:
    out = []
    for m in re.finditer(r"^[ \t]*(?:theorem|lemma)\s+(\w+)", text, re.MULTILINE):
        name = m.group(1)
        depth = 0
        i = m.end()
        end = -1
        while i < len(text) - 1:  # -1 because of `text[i+1]` lookahead on next line
            ch = text[i]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == ":" and text[i + 1] == "=" and depth == 0:
                end = i
                break
            i += 1
        if end > 0:
            out.append((m.start(), end, name))
    return out


@dataclass
class Mutation:
    module_lean: str
    theorem_name: str
    original_literal: int
    new_literal: int
    mode: str  # "single" or "global"
    mutated_text: str


def build_mutation(mut: Mutation) -> tuple[Mutation, bool]:
    fd, tmp_path = tempfile.mkstemp(
        prefix="fuzz_", suffix=".lean", dir=tempfile.gettempdir()
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
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception:
        # Treat infrastructure failures (timeout, OS errors) as build-broken,
        # so a flaky toolchain cannot manufacture spurious FREE findings.
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


def _is_skipped_literal(original: str, abs_pos: int) -> bool:
    preceding = original[max(0, abs_pos - 10) : abs_pos]
    if SKIP_CONTEXT.search(preceding):
        return True
    preceding_long = original[max(0, abs_pos - 30) : abs_pos]
    if SKIP_AFTER_INDEX.search(preceding_long):
        return True
    line_start = original.rfind("\n", 0, abs_pos) + 1
    line_prefix = original[line_start:abs_pos]
    if "--" in line_prefix:
        return True
    return False


def collect_mutations(module_lean: str) -> list[Mutation]:
    path = LEAN_ROOT / module_lean
    original = path.read_text(encoding="utf-8")
    sigs = find_signatures(original)
    muts: list[Mutation] = []
    for hstart, hend, name in sigs:
        header = original[hstart:hend]

        # ----- Mode 1: single-position mutations -----
        for lm in LITERAL.finditer(header):
            literal = int(lm.group(1))
            if literal > 99:
                continue
            abs_pos = hstart + lm.start()
            if _is_skipped_literal(original, abs_pos):
                continue
            for delta in (+1, -1):
                new_val = literal + delta
                if new_val < 0:
                    continue
                new_lit = str(new_val)
                mutated = (
                    original[:abs_pos]
                    + new_lit
                    + original[abs_pos + len(lm.group(1)) :]
                )
                muts.append(
                    Mutation(
                        module_lean, name, literal, new_val, "single", mutated
                    )
                )

        # ----- Mode 2: global-replacement mutations -----
        distinct_literals = set()
        for lm in LITERAL.finditer(header):
            n = int(lm.group(1))
            if n > 99:
                continue
            if _is_skipped_literal(original, hstart + lm.start()):
                continue
            distinct_literals.add(n)
        for literal in sorted(distinct_literals):
            for delta in (+1, -1):
                new_val = literal + delta
                if new_val < 0:
                    continue
                positions = []
                for lm in LITERAL.finditer(header):
                    if int(lm.group(1)) != literal:
                        continue
                    if _is_skipped_literal(original, hstart + lm.start()):
                        continue
                    positions.append(lm.span())
                if not positions:
                    continue
                new_header = header
                for s, e in reversed(positions):
                    new_header = new_header[:s] + str(new_val) + new_header[e:]
                mutated = original[:hstart] + new_header + original[hend:]
                muts.append(
                    Mutation(
                        module_lean, name, literal, new_val, "global", mutated
                    )
                )
    return muts


def _finding_key(mode: str, module: str, theorem: str, orig: int, new: int) -> tuple:
    return (mode, module, theorem, orig, new)


def load_allowlist() -> set[tuple]:
    if not ALLOWLIST_PATH.exists():
        print(f"  ERROR: allowlist not found at {ALLOWLIST_PATH}")
        print("         Cannot validate fuzzing results without allowlist.")
        sys.exit(1)
    data = json.loads(ALLOWLIST_PATH.read_text(encoding="utf-8"))
    return {
        _finding_key(e["mode"], e["module"], e["theorem"], e["orig"], e["new"])
        for e in data
    }


def main() -> int:
    print("=" * 72)
    print("CONSTANT FUZZING (parallel) - Tier 1 literal independence checker")
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
        print(f"  {tgt}: {len(muts)} literal perturbations to test")
        all_muts.extend(muts)

    print(f"\nDispatching {len(all_muts)} mutations across {WORKERS} workers...")
    findings: list[tuple[Mutation, bool]] = []
    with Pool(processes=min(WORKERS, max(1, len(all_muts)))) as pool:
        for i, (mut, ok) in enumerate(
            pool.imap_unordered(build_mutation, all_muts), start=1
        ):
            findings.append((mut, ok))
            tag = "[FREE]" if ok else "[load-bearing]"
            print(
                f"  ({i}/{len(all_muts)}) {tag} {mut.mode:6s} "
                f"{mut.module_lean}::{mut.theorem_name} "
                f"{mut.original_literal}->{mut.new_literal}"
            )

    free = [m for (m, ok) in findings if ok]
    free_keys = {
        _finding_key(m.mode, m.module_lean, m.theorem_name,
                     m.original_literal, m.new_literal)
        for m in free
    }
    allowlist = load_allowlist()

    # Note on counts: a single (mode, mod, thm, orig, new) key may correspond
    # to multiple list entries when the literal appears at multiple positions
    # in a theorem header (single-mode generates one Mutation per position).
    # The allowlist comparison is set-based, so the meaningful equality is
    # `unique FREE keys == allowlist size`. Both raw and unique counts are
    # printed below to make this explicit.
    print()
    print("=" * 72)
    print(f"  Total perturbations:    {len(findings)}")
    print(f"  Load-bearing:           {len(findings) - len(free)}")
    print(f"  FREE (raw findings):    {len(free)}")
    print(f"  FREE (unique keys):     {len(free_keys)}")
    print(f"  Allowlist (unique keys):{len(allowlist)}")

    unexpected = free_keys - allowlist
    missing = allowlist - free_keys

    if unexpected:
        print()
        print(f"UNEXPECTED FREE LITERALS ({len(unexpected)}, not in allowlist):")
        for mode, module, theorem, orig, new in sorted(unexpected):
            print(f"  [{mode}] {module}::{theorem}  {orig} -> {new}")

    if missing:
        print()
        print(f"MISSING FROM CURRENT RUN ({len(missing)}, in allowlist but "
              f"now load-bearing):")
        for mode, module, theorem, orig, new in sorted(missing):
            print(f"  [{mode}] {module}::{theorem}  {orig} -> {new}")

    print("=" * 72)
    if unexpected or missing:
        print("  RESULT: FAIL (allowlist drift)")
        return 1
    print(f"  RESULT: PASS ({len(free_keys)} unique FREE keys match allowlist exactly)")
    print("  Note: allowed FREE findings are intrinsic lift-norm symmetries")
    print("  (reflection a <-> q-a and saturation q >= 2a), not defects.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
