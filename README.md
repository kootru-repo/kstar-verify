# K* Verification System

[![CI](https://github.com/kootru-repo/kstar-verify/actions/workflows/verify.yml/badge.svg)](https://github.com/kootru-repo/kstar-verify/actions/workflows/verify.yml)
[![Heartbeat](https://github.com/kootru-repo/kstar-verify/actions/workflows/heartbeat.yml/badge.svg)](https://github.com/kootru-repo/kstar-verify/actions/workflows/heartbeat.yml)
[![doc-gen4](https://github.com/kootru-repo/kstar-verify/actions/workflows/docgen.yml/badge.svg)](https://kootru-repo.github.io/kstar-verify/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kootru-repo/kstar-verify/v1.0.0-submission?filepath=notebook%2Fk_star_demo.ipynb)
[![Lean](https://img.shields.io/badge/Lean4-v4.29.0--rc8-blueviolet)](lean4/lean-toolchain)
[![sorry-free](https://img.shields.io/badge/sorry-0-brightgreen)](lean4/generated/verification_manifest.json)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pending%20submission-lightgrey.svg)](SUBMISSION_CHECKLIST.html)
<!-- Replace the DOI badge above with the Zenodo SVG after the submission tag is pushed:
     [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.NNNNNNN.svg)](https://doi.org/10.5281/zenodo.NNNNNNN)
     See SUBMISSION_CHECKLIST.html "Human-gated: Zenodo DOI minting" for the update workflow. -->

This repository is the machine-executable companion to a paper on
**Krawtchouk spectral correspondence between Hamming schemes and
quantum measurement design**.  The paper argues that a specific 137-
operator Pauli set (K\*) reconstructs 4-qubit quantum states more
accurately than random Pauli sampling on noisy hardware; this repo
proves every mathematical statement in Lean 4 and replicates every
hardware fidelity from bundled IBM / Rigetti raw data.

Every mathematical claim and empirical result in the manuscript is verified
by deterministic code with no AI in the loop.

**Reviewer?** The fastest path is the Jupyter notebook
[`notebook/k_star_demo.ipynb`](notebook/k_star_demo.ipynb)
(~3 min end-to-end, all data bundled; runs every paper claim with
PASS/FAIL assertions).  Click the Binder badge above — nothing to
configure.  For the tiered
30s / 2min / 10min / (+2min optional QPU) / 30min / 60min onramp
(manifest peek, axiom table, notebook, live-QPU preflight, Lean
rebuild, Docker rebuild), see [VERIFICATION.md](VERIFICATION.md).

The notebook is the **demo-and-assert** surface: one executable file,
concrete outputs, loud asserts on drift.  The tiered `run_all.py`
script suite below is the **depth** surface: eight independent
computational backends (Lean 4, SageMath, SymPy, NumPy, cvxpy) producing
1,052 cross-checked verifications.

## Quick start (any OS, ~5 min, math-only)

Requires Python >= 3.10.  Docker and Lean 4 are optional
(see the [Full verification](#full-verification-with-hardware-data)
section below).  Raw QPU data is **already bundled** under `data/`
(~13 MB, 408 files — identical to the Zenodo archive), so all
hardware-data tiers run out of the box.

| Requirement | Version | What it enables |
|-------------|---------|-----------------|
| Python + pip | >= 3.10 | **Default quick-start: tiers 2, 3, 5** (243 math-only checks, SymPy + NumPy) |
| Docker (optional) | any | Tier 1: SageMath/FLINT exact arithmetic (44 checks) |
| Lean 4 / elan (optional) | any | Tier 0: formal proofs (11 paper theorems, 0 sorry) |

Tiers 4, 6, 7 (additional 747 hardware-data + manuscript-claim checks)
run against the repo's bundled `data/` directory — no download step.

```
git clone https://github.com/kootru-repo/kstar-verify.git
cd kstar-verify
```

Optionally create a virtual environment first
([uv](https://docs.astral.sh/uv/) recommended):
```
uv venv
source .venv/bin/activate          # Linux/macOS
.venv\Scripts\activate             # Windows
```

Then install and run the math-only tiers (243 checks, ~5 min, no data
download required):
```
pip install -r requirements.txt       # or: uv pip install -r requirements.txt
python run_all.py --registry proofs_registry.yaml --skip-tier1 --skip-lean4 --tier 2 3 5
```

This runs tier 2 (SymPy exact arithmetic, 160 checks), tier 3 (NumPy
numerical bounds, 73 checks), and tier 5 (HS->F fidelity recovery,
10 checks).  Every check prints `[PASS]` with its computed value so a
reviewer sees independent numerical evidence, not just a final "OK".

## Full verification (with hardware data)

The data-dependent tiers (4, 6, 7) replay the paper's hardware results
against the raw QPU measurement counts bundled under `data/`.  No
download step — a clone of this repo is all you need:

```
# Run all six Python tiers (2, 3, 4, 5, 6, 7) against the bundled data.
python run_all.py --registry proofs_registry.yaml --skip-tier1 --skip-lean4
```

This adds tier 4 (418 independent-verification checks against raw QPU
data), tier 6 (44 end-to-end implication chains), and tier 7 (285
manuscript-claim checks that parse every number directly from the
.tex sources).  Total wall clock: ~15 min on a modern laptop.

The Zenodo deposit (see DOI badge above) is the canonical archive for
the paper record; advanced users can point `KSTAR_DATA_DIR` at an
alternative extracted deposit to replay against that version instead.

> **Note:** `--skip-tier1` skips the SageMath/Docker tier.  If you have Docker
> installed, omit it to run all tiers.  Add `--skip-lean4` to also skip
> the Lean 4 compilation check (the sorry audit still runs).

The full system has eight tiers (0-7); tiers 0 and 1 require Lean 4 and
Docker respectively.

For math-only verification (no hardware data checks):
```
python run_all.py --registry proofs_registry.yaml --tier 2 3 5
```

> **Windows:** use `python` (not `python3`, which may be a Store stub).

---

## Running individual tiers

Follow [Quick start](#quick-start-any-os-5-min-math-only) to install, then
run tiers individually:

### Pure-math tiers (no data needed)

```
python tier2-sympy/run_all.py                       # 160 exact checks
python tier3-numpy/run_all.py                       #  73 numerical checks
python tier5-fidelity/test_fidelity_recovery.py     #  10 fidelity checks
```

Each prints PASS/FAIL per check and a summary. Any FAIL is a genuine
verification failure tied to a specific manuscript claim.

### Data-consistency tiers (QPU data bundled)

Tiers 4 and 6 verify manuscript claims against raw QPU data. The data
files are bundled in `data/` and auto-discovered by `run_all.py`.

```
python tier4-independent/run_all.py        # 418 data-consistency checks
python tier6-chains/test_implication_chain.py   #  44 end-to-end checks
python tier7-claims/run_all.py             # 285 manuscript-claim checks
```

Tier 7 parses every numerical value directly from the LaTeX manuscript,
supplemental material, and cover letter, then verifies each against its
authoritative data source (QPU JSON, simulation output, or algebraic
recomputation). Any number changed in the `.tex` without updating the
data (or vice versa) is caught.

To use an alternative data directory (e.g., from the Zenodo archive):
```
export KSTAR_DATA_DIR=/path/to/data        # Linux/macOS
set KSTAR_DATA_DIR=C:\path\to\data         # Windows (cmd)
```

### SageMath via Docker (optional)

Re-verifies Tier 1 using SageMath's FLINT/GMP backend, independent of
Python's numeric stack.

```
docker build -t kstar-verify .
docker run --rm kstar-verify
```

### Lean 4 formal proofs (optional)

```
cd lean4 && lake build && cd ..      # ~60 min first build, 0 sorry
```

#### Tier-interp: machine-checked audit of the formal statements

The Lean tier proves theorems; `tier-interp/` audits the *statements
themselves* against three failure modes that machine-checked proofs
alone cannot catch:

| Script | What it catches | Expected output |
|--------|-----------------|------------------|
| `hypothesis_ablation.py` | Decorative hypotheses (binders the proof never used). Drops each hypothesis one at a time and re-builds. | 18/18 load-bearing, 0 decorative |
| `constant_fuzzing.py` | Numeric-literal independence. Perturbs every integer constant in each Tier 1 theorem signature by +/-1 (single-position and global modes) and re-builds. Surviving perturbations are checked against `tier-interp/free_allowlist.json`, which records 69 intrinsic survivors traceable to lift-norm function symmetries (reflection a <-> q-a and saturation q >= 2a). Drift in either direction fails. | 69 FREE = allowlist, exit 0 |
| `lemma_completeness_audit.py` | Drift between manuscript labels and Lean theorem names. Walks every `\label{lem|thm|prop|cor:...}` in the .tex sources and confirms each is either mapped to an existing Lean theorem or explicitly classified as informal. Reports UNCLASSIFIED, ROTTED (manifest -> missing Lean), STALE (manifest -> dropped tex label). | 0 unclassified, 0 rotted, 0 stale |

```
cd tier-interp
python hypothesis_ablation.py        # ~10 min, parallel
python constant_fuzzing.py           # ~30 min, parallel
python lemma_completeness_audit.py   # <5 sec, no build needed
```

All three scripts auto-resolve `LEAN_ROOT` to `../lean4` relative to
the script. Override with `KSTAR_LEAN_ROOT=/path/to/lean4` if needed.
Worker count is `KSTAR_WORKERS` (default 16).

### One-command alternatives

```
# Cross-platform (Python, any OS):
python run_all.py --registry proofs_registry.yaml --tier 2 3 5

# All tiers including hardware data (bundled in data/, auto-discovered):
python run_all.py --registry proofs_registry.yaml

# Shell script (Linux/macOS/WSL only):
./validate.sh --quick       # Python tiers only
./validate.sh               # Full: Lean4 + Docker + Python
```

### Interpreting results

| Result | Meaning |
|--------|---------|
| `N/N passed` | All checks in that tier passed |
| `FAIL: <claim_id>` | A specific manuscript claim failed verification |
| `SKIP` | Tier skipped (missing data, Docker, or Lean toolchain) |

| Scope | Checks | Cumulative | Requires |
|-------|--------|------------|----------|
| Tiers 2 + 3 + 5 (math only) | 243 | 243 | Python + pip |
| + Tiers 4 + 6 (hardware data) | +462 | 705 | + QPU JSON files (bundled) |
| + Tier 7 (manuscript claims) | +285 | 990 | + LaTeX source (bundled) |
| + Tier 1 (SageMath) | +44 | 1034 | + Docker |
| + Tier 0 (Lean 4) | +18 | 1052 | + elan toolchain |

> Tiers are grouped by dependency, not numbered sequentially: tiers 2, 3, 5
> need only Python; tiers 4, 6 also need QPU data; tier 1 needs Docker;
> tier 0 needs the Lean 4 toolchain.

### What the output proves

Every check computes a value from scratch and compares it against
`proofs_registry.yaml` (the single source of truth for all manuscript
claims, SHA-256 hashed in the report). The output shows the computed
value, not just PASS/FAIL. For example:

```
  [PASS] lambda_0 = 144  (got 144)
  [PASS] N_4(5) = 137 (incl. origin)  (got 137)
  [PASS] eps_pos(W, k=2) = 11/256 = 0.042969  (registry: 11/256)
```

At the end of a `run_all.py` run, the **registry digest** prints the
key manuscript values that were verified, and the **claim coverage**
section maps every formal claim to the tiers that independently
verified it. Three report formats are generated in `results/`:

- **HTML** (`latest.html`): color-coded, collapsible report with
  executive summary, manuscript-value digest, and claim coverage
  matrix. Open in any browser.
- **Text** (`latest.txt`): full per-check details.
- **JSON** (`latest.json`): machine-readable for automated pipelines.

For a static view of which claims are covered by which tiers, see
[`coverage_matrix.md`](coverage_matrix.md).

### Proving non-vacuity (mutation testing)

To confirm the system genuinely catches errors (not just printing
PASS unconditionally), run the built-in canary mode:

```
python run_all.py --canary --registry proofs_registry.yaml
```

This injects four known-wrong values (e.g., N_4(5) = 138 instead of
137) and verifies that each mutation is detected as a FAIL. Output:

```
CANARY RESULTS: 4/4 mutations detected
  All injected errors were caught
```

### Artifact integrity

Every data file, `.tex` source, and bundled script is SHA-256 pinned in
`data/checksums.json`. This catches silent corruption, stale copies, and
copy errors between the three archives (kstar-verify, DOI, APS).

```
# CI mode: verify all 69 pinned hashes (runs in GitHub Actions)
python verify_integrity.py

# Local dev: also cross-check against DOI and APS directories
python verify_integrity.py --cross-archive

# Regenerate manifest after intentional data updates
python verify_integrity.py --generate-manifest
```

## Architecture

```
kstar-verify/
├── lean4/               Tier 0: Lean4/Mathlib4 formal proofs (0 sorry)
├── tier1-sage/          Tier 1: SageMath exact arithmetic (44 tests)
├── tier2-sympy/         Tier 2: SymPy exact rational arithmetic (160 checks)
├── tier3-numpy/         Tier 3: NumPy numerical bounds (73 checks)
├── tier4-independent/   Tier 4: Independent replication (418 checks)
├── tier5-fidelity/      Tier 5: HS-to-fidelity conversion (10 checks)
├── tier6-chains/        Tier 6: End-to-end implication chains (44 checks)
├── tier7-claims/        Tier 7: Manuscript claim verification (285 checks)
├── tier-interp/         Statement audits: ablation, fuzzing, completeness
├── notebook/            Jupyter demo (`k_star_demo.ipynb`, 74 cells, <3 min)
├── data/                QPU JSON files (408 files, ~13 MB; bundled in repo)
├── results/             Generated reports (HTML, text, JSON)
├── run_all.py              Master orchestrator
├── verify_integrity.py     Artifact integrity (SHA-256 + cross-archive)
├── registry.py             Loads claim values from proofs_registry.yaml
├── proofs_registry.yaml    Manuscript claims (single source of truth)
├── data/checksums.json     SHA-256 manifest for all data files
├── coverage_matrix.md      Every claim x every tier
├── Dockerfile              SageMath/Tier 1 container
└── requirements.txt        Python dependencies
```

| Tier | Backend | Checks | What it verifies |
|------|---------|--------|------------------|
| 0 | Lean4/Mathlib4 | 18 | Proof-theoretic endpoints (0 sorry) |
| 1 | SageMath (FLINT/GMP) | 44 | Exact arithmetic: eigenvalues, combinatorics, purity bounds |
| 2 | SymPy (Python) | 160 | Exact arithmetic: Krawtchouk, Gram, basin separation, universality |
| 3 | NumPy (Python) | 73 | Numerical: MLE convergence, Hessian, fidelity bounds, Delsarte |
| 4 | NumPy+SymPy+cvxpy | 418 | Independent: hardware data, figures, SOTA, RDM, qutrits |
| 5 | NumPy (Python) | 10 | HS-to-fidelity conversion: exact identity, near-pure recovery |
| 6 | NumPy (Python) | 44 | End-to-end: lattice-to-allocation, saturation-to-fidelity, data-to-claims |
| 7 | NumPy (Python) | 285 | Manuscript claims: parses every number from LaTeX, verifies against data |

## Independence guarantees

The verification system uses three independent computational backends:

1. **SageMath** (Tier 1): FLINT/GMP/Pari exact arithmetic -- completely
   independent of Python's numeric stack.
2. **SymPy** (Tier 2): Python exact rationals -- independent reimplementation
   of every algorithm (Krawtchouk, Gram, lattice counting).
3. **NumPy** (Tiers 3-6): IEEE 754 double-precision -- tests numerical
   predictions, MLE convergence, hardware data consistency.

Within Tier 4, seven of thirteen verification scripts import **zero
project code** -- they reimplement everything from scratch. The remaining
six import `tier4-independent/core.py` (the project's operator-selection
routines, which are the code under test) and cross-check outputs against
gold-standard libraries (cvxpy SDP).

If a reviewer suspects a SageMath-specific bug (e.g. in Krawtchouk polynomial
evaluation), Tier 2 independently reimplements every algorithm using SymPy —
a completely separate CAS backend. No single-vendor bug can produce a false pass.

## Hardware data

Tiers 4 and 6 verify hardware claims against raw QPU data bundled in
`data/` (auto-discovered by all scripts). To use an alternative data
directory (e.g., a Zenodo download), pass `--data-dir /path/to/data` or
set `KSTAR_DATA_DIR`. Tiers 2, 3, and 5 require no data files -- they
verify mathematical claims using only exact computation.

This repository bundles the complete data tree under `data/` (~13 MB,
408 files) — identical contents to the Zenodo DOI archive cited as
`\cite{zenodo}` in the manuscript, just mirrored in-repo so the full
1052-check suite runs on a cold clone with no download step.

**Note on data-consistency checks (Tier 4):** These verify reproducibility
of analysis from raw QPU JSON files — not correctness of raw measurements,
which would require re-executing circuits on hardware.

### Live-QPU replication (optional)

Anyone with an IBM Quantum account can rerun the flagship 4-qubit
W-state K* measurement on their own hardware time.  The runner,
analyzer, and acceptance bounds live under [`hardware/`](hardware/)
and are fully independent of the data-replay tiers above — no CI
cost, no bundled credentials, token read from the
`IBM_QUANTUM_TOKEN` environment variable.

```
pip install -r hardware/requirements.txt         # adds qiskit-ibm-runtime
python hardware/run_w_state.py --dry-run         # print protocol + budget
export IBM_QUANTUM_TOKEN='...'
python hardware/run_w_state.py --backend ibm_fez --n-runs 4
python hardware/analyze_results.py --reference   # recheck shipped JSON
```

The runner submits 274 measurement circuits per run (137 K* +
137 random Pauli), locally reconstructs both density matrices, and
compares F(K*), F(rand), and ΔF against the flagship reference at
[`data/w_repeat_results.json`](data/w_repeat_results.json).
See [`hardware/README.md`](hardware/README.md) for prerequisites,
budget, troubleshooting, and offline reanalysis.

## Lean4 formalization

The `lean4/` directory contains a sorry-free Lean4/Mathlib4 formalization
covering all 11 paper statements. Architecture:

- **Layer 1** (combinatorics): Pure computation via `native_decide` — verifies
  concrete n=4 instances (the paper's operational dimension). General-n properties
  are verified numerically by Tiers 2–3 (SymPy/NumPy).
- **Layer 2** (quantum linear algebra): Proved from the 5 propositional axioms listed below
- **Layer 3** (probability): Hypergeometric bounds (general n, inductive proofs)

The trust base comprises 5 propositional axioms for quantum linear
algebra (Uhlmann fidelity, Fuchs-van de Graaf complementarity, fidelity
non-negativity, Weyl eigenvalue perturbation, witness-state validity)
plus 6 opaque-type declarations for MLE/Hoeffding statistical quantities
that cannot be synthesized from the Pauli algebra. All are standard
textbook results declared in `KstarFormal/Axioms.lean`. Running
`#print axioms <theorem>` in Lean reveals exactly which axioms each result
depends on.

## What is NOT verified

| Gap | Reason |
|-----|--------|
| Natural-language proof text | Requires human mathematical judgement |
| PSD feasibility (Lemma 2(ii)) | QMA-complete -- no polynomial verification |
| M_n/N < 1 for all n | Open conjecture (explicit hypothesis in Theorem 3) |
| Trapped-ion hardware | No QPU access during submission |

### Irreducible handwritten gaps

Three proof steps in the manuscript rely on mathematical reasoning that
cannot be reduced to finite computation:

1. **Delsarte LP certificate** (Theorem 2, support-completeness): The dual
   LP feasibility argument uses Bose-Mesner idempotent structure to show no
   subset of K* operators suffices. This is a certificate-style proof over
   a continuous dual variable -- not a finite check. Lean4 axiomatizes the
   conclusion; Tiers 2-3 verify numerical consequences (eigenvalue positivity,
   Gram rank).

2. **PD Hessian implies informative-subspace identifiability** (Theorem 1):
   Positive-definiteness of the Fisher information Hessian on I_0 combined
   with global concavity of the log-likelihood implies all maximizers agree
   on informative coordinates. Non-uniqueness can only occur in flat
   (uninformative) directions, already captured by epsilon_flat. Tiers 3-4
   verify PD numerically and confirm MLE convergence empirically.

3. **Asymptotic Cramer-Rao efficiency** (Theorem 1(iii)): The claim that MLE
   achieves the Cramer-Rao bound as N_s -> infinity is a classical asymptotic
   statistics result (van der Vaart, Ch. 5). The proof is handwritten; Tier 3
   verifies the convergence rate numerically for finite N_s.
