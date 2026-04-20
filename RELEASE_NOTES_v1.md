# K* Verification Suite v1.0.0 — companion manuscript snapshot

Machine-verified certificate for the companion manuscript on Krawtchouk
spectral correspondence and measurement design.  This release is the
frozen artifact referenced by the submission.

## What this release contains

- **Lean4 formal proofs** — 22 claims kernel-checked at `v4.29.0-rc8`,
  0 `sorry`, axiom footprint audited per-claim.
- **SageMath exact-arithmetic verification** (Tier 1) — combinatorial
  identities, K\* counts, Krawtchouk eigenvalues.
- **SymPy, NumPy, independent-re-derivation tiers** (Tiers 2–4).
- **Hardware fidelity reconstructions** (Tier 5) — replay pipelines
  for IBM `ibm_fez` and Rigetti `Ankaa-3` raw data.
- **Dependency-chain verification** (Tier 6) — DAG walk over
  `depends_on` edges.
- **Manuscript-claim auditor** (Tier 7) — parses `.tex` sources and
  asserts every numerical value matches its declared data source.
- **Jupyter demo** — `notebook/k_star_demo.ipynb` runs end-to-end on
  Binder in ≈3 min; every figure and table is replayed with
  `assert` statements.

## Reviewer onramp

Three paths, time-budgeted:

| Budget | Artifact | What it proves |
| --- | --- | --- |
| 30 s | [`lean4/generated/verification_manifest.json`](https://github.com/kootru-repo/kstar-verify/blob/master/lean4/generated/verification_manifest.json) | git SHA, Lean toolchain, 0 sorry, per-theorem axiom list |
| 3 min | [Binder → `notebook/k_star_demo.ipynb`](https://mybinder.org/v2/gh/kootru-repo/kstar-verify/v1.0.0-submission?filepath=notebook%2Fk_star_demo.ipynb) | Every paper figure/table, with loud asserts |
| 45 min | `lake build` + `./validate.sh` | Independent Lean rebuild + full pipeline |

Full tiered guide: [VERIFICATION.md](https://github.com/kootru-repo/kstar-verify/blob/master/VERIFICATION.md).

## Browsable documentation

Live at **<https://kootru-repo.github.io/kstar-verify/>** — every
claim in the registry links directly to its Lean4 theorem.

## Verification snapshot at release

- Claims registered: **28**
- Lean-proved: **22** (19 `proved`, 2 `proved-from-foundational`,
  1 `proved-from-mathematical`)
- Tier-covered: **27 / 28**
- `sorry` count: **0**
- Named axioms: **4** (all declared textbook results; see the
  axiom-audit table on the landing page)

## Provenance

- All hardware data bundled with job-id manifests; no dependency on
  any author-local path.
- Every figure generated from declared data by
  [`figures/gen_figures.py`](https://github.com/kootru-repo/kstar-verify/blob/master/figures/gen_figures.py)
  (fingerprint pinned in `code_fingerprint.yaml`).
- This release is cut from commit `3796a13`.

## Integrity

```
git clone https://github.com/kootru-repo/kstar-verify
cd kstar-verify
git checkout v1.0.0-submission
./validate.sh   # writes generated/verification_manifest.json
```

A passing `validate.sh` means every claim in
[`proofs_registry.yaml`](https://github.com/kootru-repo/kstar-verify/blob/master/proofs_registry.yaml)
is discharged by its declared tier(s), the Lean library builds with
`0 sorry`, and every numeric value in the manuscript matches its data
source within declared tolerance.
