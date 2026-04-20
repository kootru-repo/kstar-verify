#!/usr/bin/env python3
"""Phase M: fresh-clone reviewer-workflow simulation.

Copies the kstar-verify tree into a temp directory (skipping
build outputs and VCS metadata), then executes the VERIFICATION.md
10-minute tier end-to-end in that copy.  Passing means a reviewer
starting from a clean checkout will see exactly what we see.

Steps in order:
  1. Copy tree to TMPDIR/kstar-verify-fresh/
  2. Verify lean-toolchain, proofs_registry.yaml, and key scripts landed.
  3. Run `python scripts/generate_all.py` inside the clone.
  4. Diff the regenerated generated/ directory against the source
     copy's generated/ directory.

Any non-trivial diff means the reviewer workflow would drift from
the submitted manifest.  Exit 0 only when diffs are zero.

Usage:
    python lean4/scripts/simulate_fresh_clone.py [--keep] [--share-cache]

    --keep          keep the temp clone for inspection (default: delete)
    --share-cache   symlink the source .lake into the clone so `lake build`
                    is incremental (seconds).  Without it, the clone is
                    fully from-scratch and `lake build` may take ~15 min.
"""
from __future__ import annotations
import argparse
import atexit
import filecmp
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

LEAN4_ROOT = Path(__file__).resolve().parent.parent
KSTAR_ROOT = LEAN4_ROOT.parent

SKIP = {
    ".lake", ".git", "__pycache__", ".venv", ".pytest_cache",
    "results", "build", "dist", "node_modules", "_scratch",
}


def copy_tree(src: Path, dst: Path) -> None:
    def ignore(_dir, names):
        return [n for n in names if n in SKIP or n.endswith(".olean")
                or n.endswith(".pyc") or n.startswith(".ipynb_checkpoints")]
    shutil.copytree(src, dst, ignore=ignore)


def run(cmd: list[str], cwd: Path, timeout: int = 1800) -> subprocess.CompletedProcess:
    print(f"    $ {' '.join(cmd)}  (in {cwd.name})")
    try:
        return subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        # Produce a CompletedProcess-shaped result so callers can inspect
        # `.returncode` uniformly.
        return subprocess.CompletedProcess(
            args=cmd, returncode=124,
            stdout="", stderr=f"timeout after {timeout}s",
        )


def check_required(clone_root: Path) -> list[str]:
    required = [
        "lean4/lean-toolchain",
        "lean4/lakefile.lean",
        "lean4/lake-manifest.json",
        "lean4/KstarFormal.lean",
        "lean4/scripts/generate_all.py",
        "lean4/scripts/check_sorry.py",
        "proofs_registry.yaml",
        "VERIFICATION.md",
        "README.md",
        "CITATION.cff",
        ".zenodo.json",
        "LICENSE",
    ]
    missing = [r for r in required if not (clone_root / r).exists()]
    return missing


def diff_generated(src_gen: Path, clone_gen: Path) -> list[str]:
    """Return list of relative paths that differ."""
    if not clone_gen.is_dir() or not src_gen.is_dir():
        return ["<generated/ directory missing>"]

    diffs: list[str] = []
    for fn in ("axiom_report.json", "sm_axiom_table.tex"):
        a, b = src_gen / fn, clone_gen / fn
        if not a.is_file() or not b.is_file():
            diffs.append(f"{fn} (one side missing)")
            continue
        if not filecmp.cmp(a, b, shallow=False):
            diffs.append(fn)
    # Manifest has generated_at + git_sha which will differ across runs.
    # Compare only the stable keys.
    a_path, b_path = src_gen / "verification_manifest.json", clone_gen / "verification_manifest.json"
    if a_path.is_file() and b_path.is_file():
        a_data = json.loads(a_path.read_text())
        b_data = json.loads(b_path.read_text())
        stable_keys = (
            "schema_version", "lean_toolchain",
            "per_theorem_axioms", "named_axioms_catalog",
            "summary",  # derived counts should be deterministic given the axiom report
            # Schema 1.2 fingerprints that DO travel with the clone:
            # only `registry_sha` (the registry file is copied).
            # `constants_sha` and `tex_sha` reference the author's private
            # manuscript-source tree, which is NOT included in the clone;
            # `tier_pass_counts` needs `results/latest.json` which the
            # clone also lacks.  Including those in the comparison would
            # turn an environment difference into a spurious drift flag.
            "registry_sha",
        )
        for k in stable_keys:
            if a_data.get(k) != b_data.get(k):
                diffs.append(f"verification_manifest.json::{k}")
    return diffs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--share-cache", action="store_true")
    args = ap.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="kstar-fresh-"))
    clone = tmp / "kstar-verify-fresh"
    print(f"[simulate] cloning to {clone}")

    def _emergency_cleanup():
        # Run on interpreter exit (Ctrl-C, normal exit, etc.) as a
        # safety net in case the main flow never reached the final
        # shutil.rmtree call.  Unlink the .lake symlink first so we
        # don't accidentally traverse into the source tree.
        # Respect --keep: if the user explicitly asked to preserve the
        # clone for inspection, do nothing here.  --keep is for
        # DELIBERATE retention; crash cleanup is for the ACCIDENTAL
        # case where the main flow didn't finish.
        if args.keep:
            return
        if not tmp.exists():
            return
        dst_lake = clone / "lean4" / ".lake"
        if dst_lake.is_symlink():
            try:
                dst_lake.unlink()
            except OSError:
                pass
        try:
            shutil.rmtree(tmp, ignore_errors=True)
        except OSError:
            pass

    atexit.register(_emergency_cleanup)

    copy_tree(KSTAR_ROOT, clone)

    missing = check_required(clone)
    if missing:
        print(f"[simulate] MISSING: {missing}")
        return 2

    clone_lean4 = clone / "lean4"

    if args.share_cache:
        # NOTE: --share-cache creates a symlink into LEAN4_ROOT/.lake inside
        # a temp directory.  On single-user dev workstations and CI runners
        # this is safe, but in a shared multi-user environment (e.g. a HPC
        # node with world-writable /tmp), another user could in principle
        # race to replace the symlink target.  No code is loaded from the
        # target (only Lean .olean files), but treat `--share-cache` as a
        # trust-local-root flag.
        src_lake = LEAN4_ROOT / ".lake"
        dst_lake = clone_lean4 / ".lake"
        if src_lake.is_dir() and not dst_lake.exists():
            print(f"[simulate] symlinking .lake from source (incremental build)")
            try:
                dst_lake.symlink_to(src_lake, target_is_directory=True)
            except OSError as e:
                print(f"  symlink failed ({e}); falling back to copy (slow)")
                shutil.copytree(src_lake, dst_lake, symlinks=True)

    print("[simulate] running lake build in clone...")
    r = run(["lake", "build"], clone_lean4)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-2000:])
        print("[simulate] FAIL: lake build failed in clone")
        return 1

    print("[simulate] running generate_all.py in clone")
    r = run([sys.executable, "scripts/generate_all.py"], clone_lean4)
    print(r.stdout[-2000:])
    if r.returncode != 0:
        print(r.stderr[-2000:])
        print("[simulate] FAIL: generate_all.py failed in clone")
        return 1

    print("[simulate] diffing generated/ against source-of-truth")
    diffs = diff_generated(
        LEAN4_ROOT / "generated",
        clone_lean4 / "generated",
    )
    if diffs:
        print("[simulate] DRIFT:")
        for d in diffs:
            print(f"  - {d}")
        return 1
    print("[simulate] OK: generated/ artifacts match byte-for-byte "
          "(stable keys only)")

    if not args.keep:
        # Unlink any .lake symlink first to avoid rmtree recursing into
        # the source tree (fatal on Windows, dangerous everywhere).
        dst_lake = clone_lean4 / ".lake"
        if dst_lake.is_symlink():
            dst_lake.unlink()
        def _force_delete(func, path):
            # Read-only files (git pack metadata) need chmod before delete.
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except OSError:
                pass
        # Python 3.12+ deprecated `onerror` in favour of `onexc`.  The two
        # callbacks have different signatures, but we only use (func, path).
        if sys.version_info >= (3, 12):
            shutil.rmtree(tmp, onexc=lambda func, path, _exc: _force_delete(func, path))
        else:
            shutil.rmtree(
                tmp,
                onerror=lambda func, path, _excinfo: _force_delete(func, path),
            )
        print(f"[simulate] cleaned up {tmp}")
    else:
        print(f"[simulate] kept clone at {clone}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
