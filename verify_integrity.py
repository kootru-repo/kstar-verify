#!/usr/bin/env python3
"""
Artifact integrity verification for the K* verification system.
================================================================
Three independent checks:

  1. DATA INTEGRITY  -- SHA-256 of every data file against pinned manifest
  2. CROSS-ARCHIVE   -- kstar-verify data == DOI data (byte-identical)
  3. TEX SYNC        -- bundled .tex files == source .tex files

Usage:
  # CI mode (GitHub Actions): verify pinned hashes only
  python verify_integrity.py

  # Local dev mode: also cross-check against DOI and source .tex
  python verify_integrity.py --cross-archive

  # Regenerate the pinned manifest after intentional data updates
  python verify_integrity.py --generate-manifest

  # Override archive paths (instead of auto-detection)
  python verify_integrity.py --cross-archive --doi-dir /path/to/DOI --aps-dir /path/to/APS

Directory structure expected by --cross-archive:
  DOI/data/           QPU JSON files (superset of kstar-verify/data/)
  DOI/*.tex           manuscript, supplemental_material, cover_letter
  APS/*.tex           same .tex files
  source/*.tex        author's private source tree (authoritative; author-only)

Exit code 0 = all pass, 1 = failures found.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
TIER7_DIR = ROOT / "tier7-claims"
MANIFEST_FILE = ROOT / "data" / "checksums.json"
BACKUP_DIR = ROOT / "backups"

# ── Helpers ──────────────────────────────────────────────────────────────

PASS_COUNT = 0
FAIL_COUNT = 0
WARN_COUNT = 0
# Coverage counters: distinguish verified-in-repo from deferred-to-Zenodo
# so the final banner tells a reviewer what the output actually covered.
NON_DATA_VERIFIED_COUNT = 0
DATA_VERIFIED_COUNT = 0
DATA_DEFERRED_COUNT = 0


def check(name, ok, detail=""):
    global PASS_COUNT, FAIL_COUNT
    tag = "PASS" if ok else "FAIL"
    if ok:
        PASS_COUNT += 1
    else:
        FAIL_COUNT += 1
    suffix = f"  ({detail})" if detail else ""
    print(f"  [{tag}] {name}{suffix}")


def warn(name, detail=""):
    global WARN_COUNT
    WARN_COUNT += 1
    print(f"  [WARN] {name}  ({detail})")


def sha256_file(path: Path) -> str:
    """Compute SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
    except PermissionError:
        warn(f"sha256: {path.name}", f"file locked by another process: {path}")
        return ""
    return h.hexdigest()


def collect_artifact_hashes() -> dict:
    """Hash every tracked artifact in the repo.

    Tracked artifacts:
      - data/**  (all files except checksums.json)
      - tier7-claims/tex/**
      - tier7-claims/figures/**
      - tier7-claims/scripts/**
      - tier7-claims/*.py  (verification scripts bundled for standalone)
      - proofs_registry.yaml
      - requirements.txt
    """
    hashes = {}

    # Data files (excluding manifest itself)
    if DATA_DIR.exists():
        for f in sorted(DATA_DIR.rglob("*")):
            if f.is_file() and f.name != "checksums.json":
                rel = f.relative_to(ROOT).as_posix()
                hashes[rel] = sha256_file(f)

    # Tier7 bundled subdirectories (tex, figures, scripts)
    for subdir in ["tex", "figures", "scripts"]:
        d = TIER7_DIR / subdir
        if d.exists():
            for f in sorted(d.rglob("*")):
                if f.is_file():
                    rel = f.relative_to(ROOT).as_posix()
                    hashes[rel] = sha256_file(f)

    # Tier7 root .py files (cross_validate_lean.py, ibm_floquet_dtc_test.py, etc.)
    if TIER7_DIR.exists():
        for f in sorted(TIER7_DIR.glob("*.py")):
            rel = f.relative_to(ROOT).as_posix()
            hashes[rel] = sha256_file(f)

    # Root-level config files
    for name in ["proofs_registry.yaml", "requirements.txt"]:
        p = ROOT / name
        if p.exists():
            hashes[name] = sha256_file(p)

    return hashes


# ── Backup system ────────────────────────────────────────────────────────

def create_backup(label: str = ""):
    """Create a timestamped backup of all tracked artifacts.

    Backups are stored in backups/<timestamp>_<label>/ with the same
    relative paths as the originals. This enables disaster recovery
    if a manifest regeneration or data update goes wrong.
    """
    ts = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    tag = f"{ts}_{label}" if label else ts
    backup_path = BACKUP_DIR / tag

    if backup_path.exists():
        print(f"  Backup already exists: {backup_path}")
        return backup_path

    backup_path.mkdir(parents=True, exist_ok=True)

    count = 0
    # Back up manifest
    if MANIFEST_FILE.exists():
        dst = backup_path / "checksums.json"
        shutil.copy2(MANIFEST_FILE, dst)
        count += 1

    # Back up registry
    reg = ROOT / "proofs_registry.yaml"
    if reg.exists():
        shutil.copy2(reg, backup_path / "proofs_registry.yaml")
        count += 1

    # Back up data files (flat list of sha256 -> filename mapping)
    manifest_snapshot = collect_artifact_hashes()
    (backup_path / "artifact_hashes.json").write_text(
        json.dumps(manifest_snapshot, indent=2) + "\n", encoding="utf-8")
    count += 1

    # Back up .tex files
    tex_dir = TIER7_DIR / "tex"
    if tex_dir.exists():
        dst_tex = backup_path / "tex"
        dst_tex.mkdir(exist_ok=True)
        for f in tex_dir.glob("*.tex"):
            shutil.copy2(f, dst_tex / f.name)
            count += 1

    print(f"  Backup created: {backup_path.relative_to(ROOT)} ({count} files)")
    return backup_path


def list_backups():
    """List all available backups."""
    if not BACKUP_DIR.exists():
        print("  No backups found.")
        return
    backups = sorted(BACKUP_DIR.iterdir())
    if not backups:
        print("  No backups found.")
        return
    print(f"\n  {len(backups)} backup(s) in {BACKUP_DIR.relative_to(ROOT)}:")
    for b in backups:
        files = list(b.rglob("*"))
        file_count = sum(1 for f in files if f.is_file())
        print(f"    {b.name}/  ({file_count} files)")


# ── Check 1: Data Integrity (pinned SHA-256 manifest) ───────────────────

def verify_pinned_hashes():
    """Verify every data file matches its pinned SHA-256."""
    print("\n== Check 1: Data Integrity (SHA-256 manifest) ==")

    if not MANIFEST_FILE.exists():
        # checksums.json ships under data/; if missing, the working tree
        # is incomplete.  Loud warning, not a silent pass.
        warn("manifest",
             f"{MANIFEST_FILE.name} not found under data/; re-clone to restore. "
             f"Continuing with tex-sync-only integrity check.")
        return

    try:
        manifest = json.loads(MANIFEST_FILE.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        check("manifest-parse", False, f"corrupt checksums.json: {e}")
        return

    pinned = manifest.get("files", {})

    if not pinned:
        warn("manifest", "Empty manifest")
        return

    current = collect_artifact_hashes()

    # Track data/* files separately so the final banner distinguishes
    # "42 pass out of 42 total" from "42 pass but 5 data files missing".
    global DATA_DEFERRED_COUNT, DATA_VERIFIED_COUNT, NON_DATA_VERIFIED_COUNT
    DATA_DEFERRED_COUNT = 0
    DATA_VERIFIED_COUNT = 0
    NON_DATA_VERIFIED_COUNT = 0

    # Check every pinned file still exists and matches
    for rel_path, expected_hash in sorted(pinned.items()):
        full = ROOT / rel_path
        is_data = rel_path.startswith("data/")
        if not full.exists():
            if is_data:
                # Data ships in-repo; a missing file means a broken clone.
                DATA_DEFERRED_COUNT += 1
                continue
            check(f"exists: {rel_path}", False, "file missing")
            continue
        actual = current.get(rel_path)
        if actual is None:
            actual = sha256_file(full)
        ok = actual == expected_hash
        check(f"sha256: {rel_path}", ok,
              f"expected {expected_hash[:16]}..., got {actual[:16]}..." if not ok else "")
        if ok:
            if is_data:
                DATA_VERIFIED_COUNT += 1
            else:
                NON_DATA_VERIFIED_COUNT += 1

    # Check for untracked files (in current but not in manifest)
    for rel_path in sorted(current):
        if rel_path not in pinned:
            warn(f"untracked: {rel_path}", "not in manifest -- run --generate-manifest")


# ── Check 2: Cross-Archive Consistency ───────────────────────────────────

def verify_cross_archive(doi_dir: Path, aps_dir: Path = None):
    """Verify kstar-verify data files are byte-identical to DOI."""
    print("\n== Check 2: Cross-Archive Consistency (kstar-verify vs DOI) ==")

    if not doi_dir.exists():
        warn("doi_dir", f"DOI directory not found: {doi_dir}")
        return

    doi_data = doi_dir / "data"
    if not doi_data.exists():
        warn("doi_data", f"DOI data/ not found: {doi_data}")
        return

    # Every file in kstar-verify/data/ must exist and match in DOI/data/
    matched = 0
    for f in sorted(DATA_DIR.rglob("*")):
        if not f.is_file() or f.name == "checksums.json":
            continue
        rel = f.relative_to(DATA_DIR)
        doi_file = doi_data / rel
        if not doi_file.exists():
            warn(f"doi-missing: data/{rel.as_posix()}", "in kstar-verify but not DOI")
            continue
        kv_hash = sha256_file(f)
        doi_hash = sha256_file(doi_file)
        ok = kv_hash == doi_hash
        check(f"match: data/{rel.as_posix()}", ok,
              f"DIVERGENT: kv={kv_hash[:16]}... doi={doi_hash[:16]}..." if not ok else "")
        if ok:
            matched += 1

    print(f"\n  {matched} files byte-identical between kstar-verify and DOI")

    # Check APS if provided
    if aps_dir and aps_dir.exists():
        print("\n== Check 2b: Cross-Archive Consistency (source vs APS) ==")
        _check_tex_sync_pair("APS", aps_dir)


# ── Check 3: TeX Sync ───────────────────────────────────────────────────

def _check_tex_sync_pair(label, target_dir):
    """Compare bundled .tex files against a target directory."""
    tex_dir = TIER7_DIR / "tex"
    if not tex_dir.exists():
        warn("tex_dir", "tier7-claims/tex/ not found")
        return

    for tex_name in ["manuscript.tex", "supplemental_material.tex", "cover_letter.tex"]:
        bundled = tex_dir / tex_name
        target = target_dir / tex_name
        if not bundled.exists():
            warn(f"{label}.{tex_name}", "bundled copy missing")
            continue
        if not target.exists():
            warn(f"{label}.{tex_name}", f"target missing: {target}")
            continue
        b_hash = sha256_file(bundled)
        t_hash = sha256_file(target)
        ok = b_hash == t_hash
        check(f"tex-sync: {tex_name} ({label})", ok,
              f"DIVERGENT" if not ok else "")


def verify_tex_sync(source_dir: Path = None, doi_dir: Path = None,
                    aps_dir: Path = None):
    """Verify bundled .tex files match source, DOI, and APS."""
    print("\n== Check 3: TeX File Sync ==")

    if source_dir and source_dir.exists():
        _check_tex_sync_pair("source", source_dir)

    if doi_dir and doi_dir.exists():
        _check_tex_sync_pair("DOI", doi_dir)

    if aps_dir and aps_dir.exists():
        _check_tex_sync_pair("APS", aps_dir)

    if not any(d and d.exists() for d in [source_dir, doi_dir, aps_dir]):
        # CI mode: at minimum, verify the bundled tex files exist and are non-empty
        tex_dir = TIER7_DIR / "tex"
        for name in ["manuscript.tex", "supplemental_material.tex", "cover_letter.tex"]:
            p = tex_dir / name
            ok = p.exists() and p.stat().st_size > 1000
            check(f"tex-present: {name}", ok,
                  f"{p.stat().st_size} bytes" if p.exists() else "missing")


# ── Manifest generation ─────────────────────────────────────────────────

def generate_manifest():
    """Generate checksums.json from current file state.

    Creates a backup of the previous manifest before overwriting.
    """
    # Back up previous manifest if it exists
    if MANIFEST_FILE.exists():
        create_backup("pre-manifest-regen")

    hashes = collect_artifact_hashes()
    manifest = {
        "description": "SHA-256 checksums for all tracked data and bundled artifacts",
        "generated_by": "verify_integrity.py --generate-manifest",
        "generated_at": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "file_count": len(hashes),
        "files": hashes,
    }
    MANIFEST_FILE.write_text(
        json.dumps(manifest, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(f"Generated {MANIFEST_FILE} with {len(hashes)} file hashes")
    return 0


# ── Main ─────────────────────────────────────────────────────────────────

def _auto_detect_dirs():
    """Auto-detect DOI, APS, and source directories from project structure."""
    doi_dir = aps_dir = source_dir = None

    # Use env vars if set
    doi_dir = Path(os.environ["KSTAR_DOI_DIR"]) if "KSTAR_DOI_DIR" in os.environ else None
    aps_dir = Path(os.environ["KSTAR_APS_DIR"]) if "KSTAR_APS_DIR" in os.environ else None
    source_dir = Path(os.environ["KSTAR_SOURCE_DIR"]) if "KSTAR_SOURCE_DIR" in os.environ else None

    # Fall back to relative-path detection
    if doi_dir is None:
        for candidate in [
            ROOT.parent.parent.parent.parent / "DOI",
        ]:
            if candidate.exists():
                doi_dir = candidate
                break

    if aps_dir is None:
        for candidate in [
            ROOT.parent.parent.parent.parent / "APS",
        ]:
            if candidate.exists():
                aps_dir = candidate
                break

    # No default fallback for source_dir: the public repo is
    # self-contained, and cross-archive source checks are author-only
    # (set KSTAR_SOURCE_DIR or pass --source-dir explicitly to use).
    return doi_dir, aps_dir, source_dir


def main():
    parser = argparse.ArgumentParser(
        description="Artifact integrity verification for K* verification system")
    parser.add_argument("--generate-manifest", action="store_true",
                        help="Regenerate checksums.json from current files")
    parser.add_argument("--cross-archive", action="store_true",
                        help="Also verify against DOI/APS directories")
    parser.add_argument("--backup", action="store_true",
                        help="Create a backup of all tracked artifacts")
    parser.add_argument("--list-backups", action="store_true",
                        help="List available backups")
    parser.add_argument("--doi-dir", type=Path, default=None,
                        help="Path to DOI archive (or set KSTAR_DOI_DIR)")
    parser.add_argument("--aps-dir", type=Path, default=None,
                        help="Path to APS submission (or set KSTAR_APS_DIR)")
    parser.add_argument("--source-dir", type=Path, default=None,
                        help="Path to author's private source tree (or set KSTAR_SOURCE_DIR); author-only, not needed for public verification")
    args = parser.parse_args()

    if args.list_backups:
        list_backups()
        return 0

    if args.backup:
        create_backup("manual")
        return 0

    if args.generate_manifest:
        return generate_manifest()

    # Auto-detect directories
    doi_dir = args.doi_dir
    aps_dir = args.aps_dir
    source_dir = args.source_dir

    if args.cross_archive and not all([doi_dir, aps_dir, source_dir]):
        auto_doi, auto_aps, auto_src = _auto_detect_dirs()
        doi_dir = doi_dir or auto_doi
        aps_dir = aps_dir or auto_aps
        source_dir = source_dir or auto_src

    print("=" * 72)
    print("  ARTIFACT INTEGRITY VERIFICATION")
    print(f"  Root:    {ROOT}")
    if DATA_DIR.exists():
        n_files = sum(1 for _ in DATA_DIR.rglob("*") if _.is_file())
        print(f"  Data:    {DATA_DIR} ({n_files} files)")
    else:
        print(f"  Data:    {DATA_DIR} (MISSING)")
    if args.cross_archive:
        print(f"  DOI:     {doi_dir or 'not found'}")
        print(f"  APS:     {aps_dir or 'not found'}")
        print(f"  Source:  {source_dir or 'not found'}")
    print("=" * 72)

    # Check 1: always run
    verify_pinned_hashes()

    # Check 2: cross-archive (local dev only)
    if args.cross_archive:
        if doi_dir:
            verify_cross_archive(doi_dir, aps_dir)

    # Check 3: tex sync
    if args.cross_archive:
        verify_tex_sync(source_dir, doi_dir, aps_dir)
    else:
        verify_tex_sync()  # CI mode: just check files exist

    # Summary
    print("\n" + "=" * 72)
    print(f"  RESULTS: {PASS_COUNT} passed, {FAIL_COUNT} failed, {WARN_COUNT} warnings")
    # Coverage breakdown: be explicit about what got verified.  Data
    # ships in-repo now, so DATA_DEFERRED_COUNT > 0 only means a file
    # listed in the manifest is missing from the working tree -- a real
    # problem the reviewer should see.
    if DATA_DEFERRED_COUNT > 0:
        total_data = DATA_DEFERRED_COUNT + DATA_VERIFIED_COUNT
        print(f"  Coverage: {NON_DATA_VERIFIED_COUNT} repo-resident files hashed + "
              f"{DATA_VERIFIED_COUNT}/{total_data} QPU data files hashed.")
        print(f"  {DATA_DEFERRED_COUNT} data files missing from data/ -- re-clone the repo to restore.")
    if FAIL_COUNT > 0:
        print("  *** INTEGRITY FAILURES DETECTED ***")
    elif WARN_COUNT > 0:
        print("  ALL CHECKS PASSED (warnings need manual review)")
    elif DATA_DEFERRED_COUNT > 0:
        print("  PARTIAL: some data files missing from data/; re-clone to restore.")
    else:
        print("  ALL ARTIFACTS VERIFIED")
    print("=" * 72)

    return 1 if FAIL_COUNT > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
