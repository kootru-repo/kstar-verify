# Independent Verification Suite -- Companion Manuscript Submission

Standalone scripts that independently verify every claim in the manuscript
using gold-standard libraries. Run `python run_all.py` for the full suite
(parallel by default; use `--sequential` for serial execution).

## Scripts

| Script | What it verifies | Method | Imports project code? |
|--------|-----------------|--------|-----------------------|
| `verify_propositions.py` | Prop 1 (Gram spectral decomposition), Prop 2 (support-completeness), K=4/K=5 eigenvalues, three-sector, dimensional rigidity, eigenvalue mass | SymPy exact rational arithmetic | **No** |
| `verify_operator_set.py` | 137 operators match lattice weight budget, basis compression (29 K* / 50 random), saturation, Rigetti operator labels | Independent lattice enumeration + greedy set cover; uses `core.py` for K* selection | Yes (`core.py`) |
| `verify_hardware_fidelities.py` | All hardware fidelities from raw expectations, F(rand), Delta F, dimensional sweep, 8q compositional | cvxpy SDP (independent MLE) cross-checked against project's `robust_mle`; F(rand) from JSON | Yes (`core.py`, `robust_mle.py`) |
| `verify_ghz_reanalysis.py` | GHZ full-state F=0.496, subspace F=0.926, DFE F=0.618, XXXX coherence, entanglement witness, missing operators | Independent subspace lstsq + DFE from Pauli coefficients; standalone fidelity/Pauli generation | Minimal (`core.py` for K* labels, `robust_mle.py` as code under test for full-state MLE) |
| `verify_allocation.py` | D(S-AR) values, allocation fractions, 13.8x/28.3x ratios, three-arm hardware F(S)/F(AR)/F(UR) | Parsed from `external_simulations_output.txt`, recomputed from S/AR/UR breakdowns | **No** |
| `verify_sota_and_stats.py` | 7-strategy SOTA ranking (Table II), adaptive failure, Wilcoxon p-values, bootstrap CIs, t-test, intra-weight decomposition | Reads cached simulation JSONs; recomputes Wilcoxon from per-trial arrays (scipy) | **No** |
| `verify_rdm_fidelities.py` | 1-RDM F(K*)=0.985, 2-RDM F(K*)=0.771, full F(K*)=0.362, K* stability vs random volatility | Standalone lstsq reconstruction + independent partial trace + Uhlmann fidelity | Minimal (`core.py` for K* label set only) |
| `verify_percentile.py` | K* percentile ~97.9%, rank ~39/1890, gap vs best < 0.005 | Exhaustive enumeration of 1,890 weight-respecting subsets (n=2) + 500 statistical draws (n=4); 32-worker parallel | Yes (`core.py`, `robust_mle.py`) |
| `verify_qutrit.py` | K*=9 for q=3, qutrit Krawtchouk eigenvalues, three-arm F(K*)>F(rand) for qutrits | Independent Gell-Mann basis, q-ary Krawtchouk polynomials, qutrit MLE, three-arm simulation | **No** (fully independent) |
| `verify_fisher_frame.py` | Fisher info equivalence, frame potential=35072, coherence=0, weight-class coverage P(all wt-1)~0.0004, ablation/progressive operator counts | Numerical Pauli construction + combinatorial formulas | **No** |
| `verify_figures.py` | All 7 manuscript figures: hardcoded data in gen_figures.py vs authoritative JSON/text sources | Cross-reference gen_figures.py arrays against data files | **No** |
| `verify_sm_extended.py` | Ququint eigenvalues (q=5), K*=q^2 shell table, eigenvalue mass progression, saturation/coverage, dimensional cross-checks, ordering sensitivity, n=5 scaling, saturate-and-fill decomposition, identity-weight pathology | Lattice enumeration + text parsing of external_simulations_output.txt + JSON cross-checks | **No** |
| `verify_remaining_claims.py` | Ablation F values (6 configs), progressive addition F values, GHZ hardware hedge scan (6 values), shot allocation arithmetic, MLE convergence behavior (structured vs random Frobenius norms) | **Cross-verification** (not independent): tests project MLE behavior under controlled conditions | Yes (`core.py`, `robust_mle.py`) — code under test |

## Dependencies

Install with: `pip install -r requirements.txt`

```
numpy
sympy
cvxpy
scipy
```

`scipy` is a transitive dependency of `robust_mle.py`.

## Data files used

All in `../../data/`:
- `hardware_results_ibm_fez_20260307_*.json` (4 single-run files)
- `w_repeat_results.json` (W-state 4-run)
- `sweep_ibm_fez_20260309_194458.json` (dimensional sweep d=2,4)
- `hardware_results_ibm_fez_20260307_223638_compositional.json` (8q)
- `oq_grouped_results_20260316.json` (Rigetti)
- `three_arm_hardware_ibm_fez_20260311_142656.json` (3-arm)
- `ghz_dfe_results.json` (GHZ reanalysis)
- `external_simulations_output.txt` (simulations 1-11)
- `sota_adaptive_n4.json` (adaptive strategy comparison)
- `saturate_fill_results.json` (intra-weight decomposition)
- `bootstrap_1000_results_20260311_153823.json` (single-run bootstrap CIs)
- `rigetti_bootstrap_results.json` (Rigetti 1000-resample bootstrap)
- `rdm_fidelity_results.json` (RDM fidelity verification output)
- `percentile_results.json` (percentile ranking verification output)
- `qutrit_verification_results.json` (qutrit extension verification output)

In `../`:
- `sota_comparison_n4.json` (7-strategy SOTA comparison)
- `significance_results_n4.json` (Wilcoxon tests, 3 reconstructors)

## Separation of concerns

- `verify_propositions.py`, `verify_allocation.py`, `verify_sota_and_stats.py`,
  `verify_qutrit.py`, `verify_fisher_frame.py`, `verify_figures.py`, and
  `verify_sm_extended.py` import **zero** project code -- pure math / text parsing /
  cached JSON / independent implementations from scratch.
- `verify_operator_set.py` imports `core.py` (code under test) but independently
  enumerates the lattice and computes the greedy set cover.
- `verify_hardware_fidelities.py` imports `core.py` and `robust_mle.py` but
  cross-checks against cvxpy SDP (independent MLE implementation).
- `verify_ghz_reanalysis.py` uses standalone fidelity and Pauli generation;
  imports `core.py` only for K* label set and `robust_mle.py` as code under
  test for full-state MLE (subspace lstsq and DFE are independent).
- `verify_rdm_fidelities.py` uses standalone lstsq reconstruction, fidelity,
  and Pauli generation; imports `core.py` only for K* label set definition.
- `verify_percentile.py` imports `core.py` and `robust_mle.py` for subset
  evaluation but independently enumerates all weight-respecting subsets.
- `verify_remaining_claims.py` is **cross-verification** (not independent):
  it tests the project's own MLE under controlled conditions (ablation,
  progressive addition, hedge scan, convergence). Uses `core.py` and
  `robust_mle.py` as code under test.
