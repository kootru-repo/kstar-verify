# Reviewer verification onramp

Tiered check-this-yourself guide for the associated manuscript
"Provable measurement design for quantum state tomography from
Hamming-scheme spectral theory." Pick the time budget you have;
each tier is self-contained and answers a specific question.

| Time budget | Tool needed | What it proves |
|-------------|-------------|----------------|
| [30 s](#30-s-inspect-the-machine-generated-manifest) | browser / text editor | Snapshot: git SHA, Lean toolchain, 0 sorry, per-theorem axiom list |
| [2 min](#2-min-read-the-auto-generated-axiom-table)  | browser / text editor | Per-claim axiom footprint auto-derived from `#print axioms` |
| [10 min](#10-min-run-the-jupyter-demo)               | Python 3.10+ (or Binder) | Every paper claim demonstrated end-to-end: combinatorial identities, hardware Table I byte-for-byte, Floquet DTC, SOTA ranking, dim-dep progression, plus 14 theorem-level numerical witnesses |
| [30 min](#30-min-rebuild-lean-from-scratch)          | [elan](https://github.com/leanprover/elan), ~4 GB disk | Independent Lean rebuild, regeneration of every machine-verified artifact |
| 60 min | Docker, ~4 GB disk | Same as 30-min tier inside a pinned container image |

---

## 30 s: inspect the machine-generated manifest

Open [`lean4/generated/verification_manifest.json`](lean4/generated/verification_manifest.json).
Keys to check:

- `git_sha`: the commit this manifest was generated from.
- `lean_toolchain`: the exact Lean version used (must match [`lean4/lean-toolchain`](lean4/lean-toolchain)).
- `lake_build.exit_code`: must be `0`.
- `sorry_audit.total_sorry`: must be `0`.
- `summary`: a classification of all `total_claims` into disjoint buckets.
  `claims_lean_core_only` (no propositional axioms beyond the three Lean core axioms),
  `claims_with_named_axioms` (use one or more of the 4 named textbook axioms),
  `claims_with_native_decide` (verified at concrete `n` by Lean's native decision
  procedure), `claims_with_opaque_decls_only` (depend only on opaque type or
  function declarations like `DensityMatrix` that are not propositional assumptions),
  and `claims_without_lean_theorem` (claims tracked by the Python verification
  suite only - e.g., equation-type entries or facts tested empirically).
  The invariant `claims_classified + claims_without_lean_theorem == total_claims`
  holds after every generator run.

If a claim you care about is in `per_theorem_axioms`, its value is the
literal output of `#print axioms <theorem>` as parsed by Lean. No human is in
the loop.

---

## 2 min: read the auto-generated axiom table

Open [`lean4/generated/sm_axiom_table.tex`](lean4/generated/sm_axiom_table.tex).
Every row maps one manuscript claim to one Lean theorem and shows the
decomposition of its axiom footprint:

- `Lean core only` = only `propext`, `Classical.choice`, `Quot.sound`. The
  gold standard — the claim is provable from nothing beyond Lean's own
  foundations.
- `\textsc{u}` = Uhlmann 1976 (`fidelity_orthogonal_zero`).
- `\textsc{fvdg}` = Fuchs and van de Graaf 1999 (`fidelity_witness_traceproj_decomp`).
- `\nativedec{N}` = `N` native_decide certificates — concrete-n numerical
  verification delegated to Lean's reflection engine. Only used on
  theorems whose statement names a concrete `n`.
- `\opaquedecl{N}` = `N` opaque-type declarations (`DensityMatrix`, `qfidelity`,
  `witnessPlus`, `witnessMinus`). These are *declarations*, not
  *propositional axioms*: Lean types, not assumptions.

Regenerate the table on a fresh clone via `python lean4/scripts/generate_all.py`.

---

## 10 min: run the Jupyter demo

Open [`notebook/k_star_demo.ipynb`](notebook/k_star_demo.ipynb) locally
or on [Binder](https://mybinder.org/v2/gh/kootru-repo/kstar-verify/v1.0.0-submission?filepath=notebook%2Fk_star_demo.ipynb)
(Binder reads [`.binder/`](.binder/) for the environment).
69 cells total (~33 markdown narrative + ~36 executable code), ~3 min
end-to-end on Binder.  All data ships in-repo under `data/`; nothing
to configure.  Every paper claim has an in-notebook witness; every
assertion fires loudly on drift.

**Sections 1-6 (structural / spectral):**
- Per-claim axiom footprint parsed from `lean4/generated/`
- Universal-n closed forms: $c_0(5,n)=1+2n$, $c_1=2n(2n-1)$, $c_2=2n(n-1)$
- Krawtchouk eigenvalue table for $n \in [4..8]$ matching Lean `native_decide`
- Phase G compositional patch coverage demonstrated at $n=8$
- $K^*=5$ saturation witnesses at every $n$
- Paper <-> Lean crosslink table for the six universal-n theorems

**Section 7 (hardware):**
- 7.1: Noiseless MLE sanity check (K\* near-unity, random-137 IC deficit)
- 7.2: 10-seed tunable-noise K\*-vs-random ensemble (clearly labelled as
  pure-sim; K\* wins 10/10 at p=0.02, $\Delta F \approx +0.14$)
- 7.3: Byte-for-byte replication of Table I's **$F(K^*) = 0.872 \pm 0.021$**
  from the W-state 4-run ibm_fez data
- 7.4: Every Table I row: Product $|+\rangle^4$ (0.996), Bell (0.518),
  GHZ raw (0.496), W single-run (0.891), all byte-for-byte
- 7.4.1: Four GHZ reconstruction variants: raw 0.496, auto-switch 0.887,
  subspace 0.926, DFE 0.618
- 7.4.2: Rigetti Ankaa-3 W-state $F(K^*) = 0.815$
- 7.4.3: Advantage magnitudes: ibm_fez W $\Delta F = +0.33$, Rigetti W
  $\Delta F = +0.248$, IBM 8-qubit compositional mean $F = 0.997$,
  Rigetti 8-qubit entangled $\Delta F = +0.534$ (5/5 patches positive)
- 7.5: Floquet DTC on ibm_kingston: $\Delta_{\mathrm{DTC}} = 0.826$,
  Pearson $r = 0.94$, 8/9 sign agreement, $C_{\ell_1} = 2.02$, negativity
  $= 0.106$, period-2 stroboscopic flip (5 hardware claims from \S sec:dtc)
- 7.6: Table \ref{tab:sota} 7-strategy ranking ($K^*$ first, +0.022 over
  D-optimal) plus 2 adaptive out-of-regime controls
- 7.7: Dim-dep $\Delta(S{-}AR)$ 13.8$\times$ drop from $n{=}2$ to $n{=}5$,
  $n{=}5$ scaling ($K^*$ vs WA-random gap narrows to +0.004), ibm_fez
  three-arm hardware sweep
- 7.8: Trust-anchor documentation for the raw counts

**Section 8 (theorem numerics, 14 subsections):**
- 8.1 Lemma 2 purity bound: 200 random mixed $\rho$, assert $\sum\lambda_i^2 \le 1$
- 8.2 Lemma 3 eigenvalue monotonicity: sweep $K \in [3..8]$, pairwise assertion
- 8.3 Lemma 4 + Theorem 1(ii) coupon-collector: `scipy.stats.hypergeom`, assert
  $P(\text{all 12 wt-1 survive}) < 5\times 10^{-4}$ and E[missing] $\approx 5.55$
- 8.4 Lemma 5 Bose-Mesner iff $q \le 3$ via lift-norm dichotomy on residues
- 8.5 fact:kfull_scaling: $K_{\mathrm{full}}(n) = n$ via Gram-rank transition
- 8.6 Lemma 1 Pauli-Bloch informative count on $|W_4\rangle$ (56 informative
  ops, $\kappa_{\mathrm{info}} = 1$, pure-direction ZZZZ excluded)
- 8.7 Lemma 6 support completeness at $K^*$: $\lambda_w > 0$ for $w=0..4$
- 8.8 Theorem 1(i) Pauli-axis adversarial witness: $\rho^\pm = (I \pm P_0)/d$
  indistinguishable on $K^*$, $F(\rho^+,\rho^-) = 0$
- 8.9 Theorem 1(iii) HS error bound (central claim): three-strategy
  comparison showing $K^*$ achieves HS$^2 = 0.406$ vs weight-$\le 2$
  oracle $0.6875 = d\varepsilon_{\mathrm{pos}}$ vs random-137 $0.5625$
- 8.9.1 Theorem 1(iii) finite-shot $O(d/N_s)$ scaling via 30-seed
  Hoeffding sweep at $N_s \in \{10^2..10^5\}$ (measured/predicted ratio 0.996)
- 8.10 Theorem 2: $c_1$ jump $8 \to 56$, Delsarte zero $c_4(K{=}3) = 0$,
  $\kappa(G) = 256/64 = 4$
- 8.11 Theorem 3 asymptotic: $M_n/N_n$ ratios for $n{=}7,8,9$ match registry
- 8.12 Corollary 1 approximate locality: 41 high-weight nonzero ops on
  $|W_4\rangle$, $\sum x_P^2 = 11$, $\varepsilon_{\mathrm{tail}} = 11/256$
- 8.13 Corollary 2 operator lower bound: $\sum \binom{n}{w}(q^2{-}1)^w = 66$,
  $K^*$ weight-$\le 2$ count $= 67$ (margin $= 1$ for identity)

Nothing to configure — the repo's `data/` directory is auto-discovered
and all 408 hardware-data files (~13 MB, identical to the Zenodo
archive) are bundled.  Advanced users
replaying against an alternative deposit can set `KSTAR_DATA_DIR` to
override the auto-discovered path.

Every cell exits deterministically with PASS/FAIL assertions.  Any
drift between the registry and the computed value fires a loud
`AssertionError` identifying the claim.

---

## +2 min (optional): live-QPU preflight

**If you have a free IBM Quantum account, you can re-run the W-state
K\* protocol on real hardware yourself.**  The bundled Table I
replication above answers *"do the claimed numbers come out of the
raw counts?"*; this step answers the stronger question *"does the
same protocol still produce a K\* advantage on hardware today?"*

```
pip install -r hardware/requirements.txt
export IBM_QUANTUM_TOKEN='your-token-here'

python hardware/run_w_state.py --verify-equivalence   # 0 QPU
python hardware/run_w_state.py --dry-run              # 0 QPU
python hardware/run_w_state.py --preflight            # ~2 s QPU
```

Each step gates the next, so a broken token / stale toolchain / or
K\* operator drift fails in seconds rather than minutes.  Preflight
submits exactly 2 measurement circuits (checks `<IIIZ>` and `<IIXI>`
on the 4-qubit W state) and verifies they fall within the expected
noise tolerance.  A `[PASS]` confirms the token, backend routing,
`SamplerV2` round-trip, counts extraction, and basis-rotation are
all correct.

For the full flagship replication (4 runs × 137 K\* + 137 random
Paulis at 1000 shots, ~2 min of QPU time, well under the 10-min
free-tier monthly budget):

```
python hardware/run_w_state.py --backend ibm_fez --n-runs 4
python hardware/analyze_results.py hardware/results/w_repeat_results_*.json
```

The emitted JSON matches the shipped reference schema, including
`provenance.raw_counts` keyed by job_id and
`provenance.calibration` capturing backend state at submit time.
A 5-minute QPU-time safety gate (`--yes-over-5min` to override)
prevents accidental budget blow-outs.  See
[`hardware/README.md`](hardware/README.md) for full details and
troubleshooting.

---

## 30 min: rebuild Lean from scratch

```
git clone https://github.com/kootru-repo/kstar-verify.git
cd kstar-verify/lean4
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y
lake exe cache get           # Mathlib binary cache (~3 GB; saves ~40 min of rebuild)
lake build                   # ~5 min on an M1; ~15 min on a CI runner
pip install pyyaml           # required by the generator scripts
python scripts/generate_all.py
```

Expected output:

- `lake build` exits 0 with 1272 jobs.
- `python scripts/check_sorry.py` reports 0 sorry across all files.
- `python scripts/generate_all.py` exits 0 through all five steps
  (dump_axioms, manifest, sm_table, sync_registry --check, link-check).
- `lean4/generated/verification_manifest.json` reflects your local
  `git_sha` and matches the submission's `per_theorem_axioms` byte-for-byte.

Any divergence from the submitted manifest points to either
(a) you are on a different commit, or (b) the submission was
tampered with; both are falsifiable by `git log`.

---

## Spot-checking a specific claim

To verify a specific claim referenced in the manuscript:

1. Find its Lean theorem name in [`proofs_registry.yaml`](proofs_registry.yaml)
   (field `lean4_theorem_universal` for universal-n, else `lean4_theorem`).
2. Open a Lean REPL: `cd lean4 && lake env lean --stdin`, then
   `import KstarFormal` and `#print axioms <theorem_name>`.
3. Cross-reference with [`lean4/generated/axiom_report.json`](lean4/generated/axiom_report.json).

---

## Reporting discrepancies

If any artifact in this submission fails to match its machine-derived
form, please file the discrepancy at the GitHub issue tracker of
`kootru-repo/kstar-verify`, including the `git_sha` from your
local manifest and the divergent output.
