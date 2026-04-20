# Machine-generated verification artifacts

Files in this directory are **auto-generated** by the Phase H
scripts in [../scripts/](../scripts/). Do not edit by hand; regenerate
after any Lean source change via:

```
python scripts/generate_all.py
```

## Files

| File | Producer | Purpose |
|------|----------|---------|
| `axiom_report.json` | `scripts/dump_axioms.py` | Per-theorem axiom list parsed from `#print axioms`. Mapping: `theorem_name -> [axiom1, axiom2, ...]`. |
| `axiom_report.txt` | `scripts/dump_axioms.py` | Raw `#print axioms` transcript. Human-auditable. |
| `verification_manifest.json` | `scripts/generate_verification_manifest.py` | Single JSON snapshot of the full Lean4 verification state: git SHA, Lean toolchain, build status, sorry count, claim-to-theorem-to-axiom map, summary statistics. Archive this with any submission. |
| `sm_axiom_table.tex` | `scripts/generate_sm_axiom_table.py` | LaTeX tabular, `\input{}`-ready, showing the axiom footprint for every claim with a registered Lean theorem. Goes into supplemental_material.tex. |

## How a reviewer uses these

1. `verification_manifest.json` tells you exactly which git commit and
   Lean toolchain produced the artifacts.
2. `axiom_report.txt` is the raw, reproducible output of the Lean
   elaborator. Re-run `python scripts/dump_axioms.py` inside a matching
   toolchain to regenerate it and diff against this copy.
3. `sm_axiom_table.tex` is what appears in the manuscript, so a
   reviewer can check that paper claims match machine evidence.

## Regeneration

Requires the Lean4 toolchain from `../lean-toolchain` on PATH (via
[`elan`](https://github.com/leanprover/elan)) and Python 3.10+ with
`pyyaml`.

```
cd lean4
python scripts/generate_all.py
```

Each script is also idempotent on its own.
