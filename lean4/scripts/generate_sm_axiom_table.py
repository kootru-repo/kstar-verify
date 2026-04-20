#!/usr/bin/env python3
"""Phase H: auto-generate the SM axiom-footprint table from the
machine-verified axiom report.

Emits a LaTeX `\\begin{tabular}` snippet to be `\\input{...}`'ed by
supplemental_material.tex under the sm:lean-universal subsection.
Columns:
  Claim ID | Lean theorem (universal) | Axiom footprint

Axiom footprint classification:
  * Lean core only       -- {propext, Classical.choice, Quot.sound}
  * + named axioms       -- one or more of the 4 named propositional axioms
  * + native_decide      -- concrete-n certificates (should appear only on
                            theorems that name a concrete n in their statement)

Inputs:
  proofs_registry.yaml
  lean4/generated/axiom_report.json   (run dump_axioms.py first)

Output:
  lean4/generated/sm_axiom_table.tex
"""
from __future__ import annotations
import json
import os
import sys
import tempfile
from pathlib import Path


def _atomic_write_text(path: Path, payload: str) -> None:
    """Atomic write via same-dir temp file + os.replace."""
    path.parent.mkdir(exist_ok=True, parents=True)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
            f.write(payload)
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise

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
# Manuscript-side mirror for \input{generated/sm_axiom_table.tex},
# relative to supplemental_material.tex.  Author-only: in a public
# clone the source tree won't exist and the copy is skipped.
SUBMISSION_ROOT = KSTAR_ROOT.parent  # submission/submission
MANUSCRIPT_GENERATED = Path(
    os.environ.get("KSTAR_MANUSCRIPT_ROOT", SUBMISSION_ROOT / "manuscript-source")
) / "generated"

LEAN_CORE = {"propext", "Classical.choice", "Quot.sound"}

NAMED = {
    "witness_states_are_valid": r"\textsc{wit}",
    "fidelity_orthogonal_zero": r"\textsc{u}",
    "fidelity_witness_traceproj_decomp": r"\textsc{fvdg}",
    "fidelity_nonneg": r"\textsc{nn}",
    "weyl_eigenvalue_bound": r"\textsc{w}",
}


def latex_escape(s: str) -> str:
    """Escape LaTeX special characters in a plain-text identifier.

    Order of substitutions matters: braces have to be escaped before
    any replacement that emits braces as part of its replacement
    text (e.g. `\\textbackslash{}`), otherwise the emitted braces
    get double-escaped.
    """
    if s is None:
        return ""
    # Braces first (protects against double-escape of \\textbackslash{}).
    out = s.replace("{", r"\{").replace("}", r"\}")
    out = out.replace("\\", r"\textbackslash{}")
    out = out.replace("_", r"\_")
    out = out.replace("&", r"\&")
    out = out.replace("%", r"\%")
    out = out.replace("$", r"\$")
    out = out.replace("#", r"\#")
    return out


def classify(axioms: list[str]) -> str:
    named_hits = [NAMED[a] for a in axioms if a in NAMED]
    has_native = any("._native.native_decide.ax_" in a for a in axioms)
    other = [a for a in axioms if a not in LEAN_CORE and a not in NAMED and "._native.native_decide.ax_" not in a]

    badges: list[str] = []
    if named_hits:
        badges.append(", ".join(sorted(set(named_hits))))
    if has_native:
        # Count native_decide certs to give a quick impression of witness count.
        ct = sum(1 for a in axioms if "._native.native_decide.ax_" in a)
        badges.append(f"\\nativedec{{{ct}}}")
    if other:
        # Non-core, non-named, non-native: usually opaque decls (DensityMatrix etc.).
        # Show count; full list is in the JSON.
        badges.append(f"\\opaquedecl{{{len(other)}}}")

    if not badges:
        return r"Lean core only"
    return "Lean core + " + "; ".join(badges)


def main() -> int:
    if not AXIOM_REPORT.is_file():
        print(f"ERROR: {AXIOM_REPORT} missing. Run dump_axioms.py first.")
        return 2
    registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    axiom_report = json.loads(AXIOM_REPORT.read_text(encoding="utf-8"))

    rows: list[tuple[str, str, str]] = []
    for stmt in registry.get("statements", []):
        cid = stmt.get("id")
        if not cid:
            continue  # Skip malformed registry entries.
        # Prefer the universal theorem when present; fall back to the
        # n=4 Lean theorem otherwise.
        tname = stmt.get("lean4_theorem_universal") or stmt.get("lean4_theorem")
        if not tname or tname not in axiom_report:
            continue
        rows.append((cid, tname, classify(axiom_report[tname])))

    # Only emit legend entries for textbook axiom codes that ACTUALLY
    # appear in the generated rows.  Unused codes (e.g., `wit`, `nn`,
    # `w` when only `u` and `fvdg` were cited) previously leaked into
    # the legend and made the table look stale to a fresh reviewer.
    AXIOM_LEGEND_ENTRIES = {
        r"\textsc{u}": r"\textsc{u}=Uhlmann~1976",
        r"\textsc{fvdg}": r"\textsc{fvdg}=Fuchs--van de Graaf~1999",
        r"\textsc{wit}": r"\textsc{wit}=Nielsen--Chuang witness validity",
        r"\textsc{nn}": r"\textsc{nn}=fidelity non-negativity",
        r"\textsc{w}": r"\textsc{w}=Weyl eigenvalue bound",
    }
    used_codes = sorted({
        code
        for code in AXIOM_LEGEND_ENTRIES
        if any(code in footprint for _, _, footprint in rows)
    })
    legend_body = "; ".join(AXIOM_LEGEND_ENTRIES[c] for c in used_codes) if used_codes else ""
    legend_line = (
        rf"% Legend: {legend_body}."
        if legend_body
        else r"% Legend: (no named textbook axioms cited in this build)"
    )

    lines: list[str] = [
        "% Auto-generated by lean4/scripts/generate_sm_axiom_table.py",
        "% DO NOT EDIT. Regenerate after any Lean change.",
        legend_line,
        r"% \nativedec{N}=N concrete-n native\_decide certs.",
        r"% \opaquedecl{N}=N opaque declarations (DensityMatrix/qfidelity-family).",
        r"% \providecommand lets this snippet be \input'd twice (e.g.,"
        r" main + appendix) without 'already defined' errors.",
        r"\providecommand{\nativedec}[1]{\textsc{nd}\textsubscript{#1}}",
        r"\providecommand{\opaquedecl}[1]{\textsc{opq}\textsubscript{#1}}",
        r"\begin{tabular}{@{}lll@{}}",
        r"\toprule",
        r"Claim ID & Lean theorem & Axiom footprint \\",
        r"\midrule",
    ]
    for cid, tname, footprint in rows:
        lines.append(
            rf"\texttt{{{latex_escape(cid)}}} & \texttt{{{latex_escape(tname)}}} & {footprint} \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]

    out = GENERATED / "sm_axiom_table.tex"
    payload = "\n".join(lines) + "\n"
    _atomic_write_text(out, payload)
    print(f"[sm_axiom_table] wrote {out}  ({len(rows)} rows)")

    # Mirror to the manuscript-side generated/ directory so SM can
    # `\input{generated/sm_axiom_table.tex}` without crossing repo roots.
    mirror = MANUSCRIPT_GENERATED / "sm_axiom_table.tex"
    if MANUSCRIPT_GENERATED.parent.is_dir():
        try:
            _atomic_write_text(mirror, payload)
            print(f"[sm_axiom_table] mirrored to {mirror}")
        except (OSError, NotADirectoryError) as e:
            # Do NOT silently succeed: a stale mirror would cause the
            # manuscript PDF to embed outdated axiom data.
            print(
                f"ERROR: primary write to {out} succeeded but mirror "
                f"to {mirror} failed: {e}"
            )
            return 1
    else:
        print(
            f"[sm_axiom_table] manuscript tree absent at "
            f"{MANUSCRIPT_GENERATED.parent}; skipping mirror."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
