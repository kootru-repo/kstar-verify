#!/usr/bin/env python3
"""Phase H: machine-generated verification manifest.

Produces a single JSON document that snapshots the complete Lean4
verification state for the associated manuscript submission.  Designed to be
uploaded as a CI artifact and archived at the Zenodo DOI so a reviewer
can verify the claim-to-theorem-to-axiom-chain without running Lean.

Inputs:
  proofs_registry.yaml                         -- claim <-> theorem map
  lean4/generated/axiom_report.json            -- per-theorem axiom list
                                                  (run dump_axioms.py first)
  lean4/lean-toolchain                         -- pinned toolchain
  lake build                                   -- build status
  scripts/check_sorry.py (Layer 1 mode)        -- sorry audit

Output:
  lean4/generated/verification_manifest.json
"""
from __future__ import annotations
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    import yaml
except ModuleNotFoundError:
    sys.stderr.write(
        "ERROR: PyYAML required.  Install: pip install pyyaml>=6.0\n"
    )
    sys.exit(2)

LEAN4_ROOT = Path(__file__).resolve().parent.parent
KSTAR_ROOT = LEAN4_ROOT.parent
REGISTRY = KSTAR_ROOT / "proofs_registry.yaml"
GENERATED = LEAN4_ROOT / "generated"
AXIOM_REPORT = GENERATED / "axiom_report.json"
TOOLCHAIN = LEAN4_ROOT / "lean-toolchain"
# Schema 1.2 fields (G9 of the CI/CD gap inventory).
# MANUSCRIPT_ROOT points at the author's private source tree when one
# is available (set KSTAR_MANUSCRIPT_ROOT to override).  In a public
# clone this directory won't exist and tex_sha / constants_sha fields
# resolve to null (handled gracefully by downstream consumers).
SUBMISSION_ROOT = KSTAR_ROOT.parent  # submission/submission
MANUSCRIPT_ROOT = Path(
    os.environ.get("KSTAR_MANUSCRIPT_ROOT", SUBMISSION_ROOT / "manuscript-source")
)
CONSTANTS_JSON = (
    MANUSCRIPT_ROOT / "scripts" / "independent-verification" / "canonical_constants.json"
)
# Intentional list: only the two source-of-claims .tex files.
# cover_letter.tex and figure-only .tex files are excluded.  Keep this
# list in sync with the SUBMISSION_CHECKLIST "Manuscript" section; adding
# a third file here changes the claim-set fingerprint (tex_sha).
TEX_FILES = [
    MANUSCRIPT_ROOT / "manuscript.tex",
    MANUSCRIPT_ROOT / "supplemental_material.tex",
]
LATEST_RESULTS = KSTAR_ROOT / "results" / "latest.json"

import tempfile


def _atomic_write_text(path: Path, payload: str, encoding: str = "utf-8") -> None:
    """Write `payload` to `path` atomically.

    Writes to a temp file in the same directory, then `os.replace`s
    onto the target.  This guarantees that a reader either sees the
    full old content or the full new content, never a truncated
    partial write (even if the process is Ctrl-C'd mid-write).
    """
    path.parent.mkdir(exist_ok=True, parents=True)
    # Same-directory temp so `os.replace` is a rename on the same filesystem.
    fd, tmp_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="\n") as f:
            f.write(payload)
        os.replace(tmp_name, path)
    except Exception:
        # Clean up the temp file if the replace didn't happen.
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


LEAN_CORE_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}

NAMED_PROPOSITIONAL_AXIOMS = {
    "witness_states_are_valid": "Nielsen & Chuang (Sec. 2.4)",
    "fidelity_orthogonal_zero": "Uhlmann 1976",
    "fidelity_witness_traceproj_decomp": "Fuchs & van de Graaf 1999",
    "fidelity_nonneg": "definition",
    "weyl_eigenvalue_bound": "Bhatia, Matrix Analysis, Thm III.2.1",
}


def git_sha() -> str:
    for root in (KSTAR_ROOT, KSTAR_ROOT.parent, KSTAR_ROOT.parent.parent):
        try:
            sha = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=root, text=True, stderr=subprocess.DEVNULL
            ).strip()
            if sha:
                return sha
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    return "local-dev"


def toolchain_version() -> str:
    try:
        return TOOLCHAIN.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        return "unknown"


def _sha256_file(path: Path) -> str | None:
    """Return SHA-256 hex digest of a file, or None if missing."""
    import hashlib
    if not path.is_file():
        return None
    h = hashlib.sha256()
    # Read in chunks so we don't page-in large .tex / JSON files.
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def compute_source_hashes() -> dict:
    """Fingerprint the files that together determine the claim set.

    Every entry is either a 64-char SHA-256 hex digest or null if the
    file doesn't exist on this machine (e.g., standalone Lean-only
    clones don't have the manuscript tree).
    """
    return {
        "registry_sha": _sha256_file(REGISTRY),
        "constants_sha": _sha256_file(CONSTANTS_JSON),
        "tex_sha": {p.name: _sha256_file(p) for p in TEX_FILES},
    }


def read_tier_pass_counts() -> dict:
    """Return per-tier pass counts from results/latest.json if available.

    Missing file, absent tiers block, or tier entries without
    `checks_passed` all resolve to an empty dict -- the manifest is
    still usable, just without the Python-tier fingerprint.
    """
    if not LATEST_RESULTS.is_file():
        return {}
    try:
        data = json.loads(LATEST_RESULTS.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    out: dict[str, int] = {}
    for tier_id, tier in (data.get("tiers") or {}).items():
        passed = tier.get("checks_passed")
        if isinstance(passed, int):
            out[tier_id] = passed
    return out


def run_lake_build() -> dict:
    """Return {exit_code, jobs}. Does NOT rebuild if artifacts are fresh."""
    try:
        proc = subprocess.run(
            ["lake", "build"], cwd=LEAN4_ROOT, capture_output=True,
            text=True, timeout=900,
        )
        out = (proc.stdout or "") + (proc.stderr or "")
        jobs = 0
        # Lake's "Build completed successfully (N jobs)" line, when present.
        for line in out.splitlines():
            if "Build completed successfully" in line:
                try:
                    jobs = int(line.split("(")[1].split()[0])
                except (IndexError, ValueError):
                    pass
        return {"exit_code": proc.returncode, "jobs": jobs}
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return {"exit_code": -1, "jobs": 0, "error": str(e)}


def run_sorry_audit() -> dict:
    """Parse check_sorry.py's Layer 1 / Layer 2-3 counts.

    Defensive: both lines must be present.  If either is absent we
    report the parse failure rather than silently returning zero
    (which would mask a missing audit).
    """
    try:
        proc = subprocess.run(
            [sys.executable, "scripts/check_sorry.py", "--layer1"],
            cwd=LEAN4_ROOT, capture_output=True, text=True, timeout=120,
        )
        layer1: int | None = None
        layer23: int | None = None
        for line in (proc.stdout or "").splitlines():
            line = line.strip()
            if line.startswith("Layer 1 sorry:"):
                layer1 = int(line.split(":", 1)[1].strip())
            elif line.startswith("Layer 2/3 sorry:"):
                layer23 = int(line.split(":", 1)[1].strip())
        if layer1 is None or layer23 is None:
            return {
                "exit_code": proc.returncode,
                "error": "check_sorry.py output format changed; "
                         "could not parse Layer 1 / Layer 2-3 counts",
                "stdout_tail": (proc.stdout or "")[-500:],
            }
        return {
            "exit_code": proc.returncode,
            "layer1_sorry": layer1,
            "layer23_sorry": layer23,
            "total_sorry": layer1 + layer23,
        }
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return {"exit_code": -1, "error": str(e)}


def classify_axioms(axioms: list[str]) -> dict:
    """Partition axioms into core / named / native_decide / other."""
    core = [a for a in axioms if a in LEAN_CORE_AXIOMS]
    named = [a for a in axioms if a in NAMED_PROPOSITIONAL_AXIOMS]
    native = [a for a in axioms if "._native.native_decide.ax_" in a]
    known = set(core) | set(named) | set(native)
    other = [a for a in axioms if a not in known]
    return {
        "core": core,
        "named_propositional": named,
        "native_decide_certificates": native,
        "other_declarations": other,
    }


def build_claim_map(registry: dict, axiom_report: dict) -> list[dict]:
    entries: list[dict] = []
    for stmt in registry.get("statements", []):
        cid = stmt.get("id")
        if not cid:
            continue
        entry: dict = {
            "claim_id": cid,
            "label": stmt.get("label", ""),
            "lean4_theorem": stmt.get("lean4_theorem"),
            "lean4_theorem_universal": stmt.get("lean4_theorem_universal"),
        }
        for key in ("lean4_theorem", "lean4_theorem_universal"):
            tname = stmt.get(key)
            if tname and tname in axiom_report:
                entry[f"{key}_axioms"] = axiom_report[tname]
                entry[f"{key}_classification"] = classify_axioms(axiom_report[tname])
        entries.append(entry)
    return entries


def main() -> int:
    if not AXIOM_REPORT.is_file():
        print(f"ERROR: {AXIOM_REPORT} missing. Run dump_axioms.py first.")
        return 2
    try:
        registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        print(f"ERROR: {REGISTRY} is not valid YAML: {e}")
        return 2
    try:
        axiom_report = json.loads(AXIOM_REPORT.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        print(
            f"ERROR: {AXIOM_REPORT} is not valid JSON: {e}\n"
            "Re-run `python scripts/dump_axioms.py` to regenerate."
        )
        return 2

    source_hashes = compute_source_hashes()
    manifest = {
        # Schema 1.3: identical to 1.2 minus the private-IP pin field
        # that was removed when the verification surface was made
        # self-contained for submission.  Every remaining field is
        # reproducible from the public tree.
        "schema_version": "1.3",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_sha": git_sha(),
        "lean_toolchain": toolchain_version(),
        "registry_sha": source_hashes["registry_sha"],
        "constants_sha": source_hashes["constants_sha"],
        "tex_sha": source_hashes["tex_sha"],
        "tier_pass_counts": read_tier_pass_counts(),
        "lake_build": run_lake_build(),
        "sorry_audit": run_sorry_audit(),
        "named_axioms_catalog": NAMED_PROPOSITIONAL_AXIOMS,
        "per_theorem_axioms": axiom_report,
        "claim_map": build_claim_map(registry, axiom_report),
    }

    # Aggregate cross-cut: how many claims prove with "Lean core only"?
    # Opaque declarations (DensityMatrix, qfidelity-family) are type/function
    # constants, not propositional axioms, so they are lumped with "core only"
    # for trust-base accounting unless the theorem also carries named axioms
    # or native_decide certs.
    core_only = 0
    core_plus_named = 0
    with_native_decide = 0
    with_opaque_only = 0
    no_lean_theorem = 0  # claims the paper makes but that have no Lean entry
    for entry in manifest["claim_map"]:
        cls = entry.get("lean4_theorem_universal_classification") or entry.get(
            "lean4_theorem_classification"
        )
        if not cls:
            # No Lean theorem at all (e.g., equation-type entry, or a fact
            # tracked only by the Python verification suite).
            no_lean_theorem += 1
            continue
        if cls["named_propositional"]:
            core_plus_named += 1
        elif cls["native_decide_certificates"]:
            with_native_decide += 1
        elif cls["other_declarations"]:
            with_opaque_only += 1
        else:
            core_only += 1
    total = len(manifest["claim_map"])
    classified = core_only + core_plus_named + with_native_decide + with_opaque_only
    manifest["summary"] = {
        "total_claims": total,
        "claims_lean_core_only": core_only,
        "claims_with_named_axioms": core_plus_named,
        "claims_with_native_decide": with_native_decide,
        "claims_with_opaque_decls_only": with_opaque_only,
        "claims_without_lean_theorem": no_lean_theorem,
        "claims_classified": classified,
    }
    # Invariant: every claim falls into exactly one bucket.
    assert classified + no_lean_theorem == total, (
        f"claim classification lost entries: {classified} + {no_lean_theorem} != {total}"
    )

    GENERATED.mkdir(exist_ok=True, parents=True)
    outpath = GENERATED / "verification_manifest.json"
    _atomic_write_text(outpath, json.dumps(manifest, indent=2) + "\n")
    print(f"[manifest] wrote {outpath}")
    print(f"[manifest] git: {manifest['git_sha'][:10]}  toolchain: {manifest['lean_toolchain']}")
    print(
        f"[manifest] build exit: {manifest['lake_build']['exit_code']}  "
        f"jobs: {manifest['lake_build']['jobs']}"
    )
    print(
        f"[manifest] sorry: layer1={manifest['sorry_audit'].get('layer1_sorry')}  "
        f"layer2/3={manifest['sorry_audit'].get('layer23_sorry')}"
    )
    s = manifest['summary']
    print(
        f"[manifest] claims: {s['total_claims']}  "
        f"core-only: {s['claims_lean_core_only']}  "
        f"named-ax: {s['claims_with_named_axioms']}  "
        f"native-decide: {s['claims_with_native_decide']}  "
        f"opaque-only: {s['claims_with_opaque_decls_only']}  "
        f"no-lean: {s['claims_without_lean_theorem']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
