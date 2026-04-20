Supplemental Material for:
"Krawtchouk-structured Pauli measurement selection outperforms
 random sampling for quantum state tomography"

Author: Andrew Craton, Kootru Labs (acraton@kootru.com)
Journal: K* Verification (submitted)

============================================================
HARDWARE DATA FILES
============================================================

IBM experiments performed on ibm_fez (156-qubit Eagle r3),
March 7, 9, 10, and 11, 2026.
Rigetti experiments performed on Ankaa-3 (84-qubit) via Open Quantum,
March 15-16 and 21-22, 2026.

File                                        State         Qubits  Reference
----                                        -----         ------  ---------
hardware_results_ibm_fez_20260307_214441    Product |+>^4  4      Table I
hardware_results_ibm_fez_20260307_214922    GHZ            4      Table I
hardware_results_ibm_fez_20260307_220545    Bell pair      2      Table I
hardware_results_ibm_fez_20260307_220810    W state        4      Table I (single run)
hardware_results_ibm_fez_20260307_223638_compositional
                                            8q composite   8      Sec. V C
sweep_ibm_fez_20260309_194458_enriched      W states n=2,4 2,4    Table I (dagger rows)
w_repeat_results_enriched                   W state x4     4      Table I (ddagger)
three_arm_hardware_ibm_fez_..._enriched     3-arm W test   4      Table III
oq_grouped_results_20260316                 Rigetti W n=4  4      Table I (Ankaa-3 rows)
8qubit-w-rigetti/                           Rigetti 8q     8      Sec. V C

============================================================
K* MEASUREMENT SET (machine-readable)
============================================================

kstar_operators_n4.json    All 137 K*-selected Pauli operators for n=4
                           Weight distribution: {0:1, 1:12, 2:54, 3:54, 4:16}
                           Deterministic: select_kstar_paulis(4) in scripts/core.py

kstar_bases_29.json        The 29 tensor-product measurement bases covering
                           all 137 operators, with per-basis operator coverage map.
                           Deterministic: greedy set-cover of operator set.

============================================================
PROVENANCE AND ANTI-SPOOFING
============================================================

provenance_manifest.json   Complete data provenance manifest (platforms,
                           job IDs, calibration, anti-spoofing evidence)
ibm_qubit_layout.json      Physical qubit maps for all IBM circuits
provenance_ibm_andrew_craton.json   IBM job metadata (andrew.craton account)
provenance_ibm_quantum_ibm_01.json  IBM job metadata (quantum_ibm_02 account)
oq_all_jobs.json           Full OQ scheduler API dump (471 jobs)

Every hardware result file contains a "provenance" block with:
  - job_ids: IBM/OQ job identifiers (retrievable with account access)
  - raw_counts: per-job bitstring measurement histograms (1000 shots each)
  - calibration / calibration_snapshot: backend gate errors at time of run
  - physical_qubits: which physical qubits were used
  - transpilation: gate set, optimization level
  - random_basis_selection: seed used for random operator selection

============================================================
ANALYSIS FILES
============================================================

hyperparameter_sweep_results.json   MLE hyperparam sweep      SM Sec. XII E
rigetti_bootstrap_results.json      Rigetti bootstrap CIs     Table I footnote
bootstrap_1000_results_20260311_153823.json
                                    IBM W-state bootstrap (1000 resamples)
ghz_dfe_results.json                GHZ reanalysis: full Hradil MLE (0.888),
                                    ratio-form MLE (0.496), subspace (0.926),
                                    DFE lower bound (0.618)         Sec. V D
mle_convergence_comparison_results.json
                                    Full vs ratio Hradil form comparison
                                    (iteration counts, fidelities)  SM Sec. VII D
sota_adaptive_n4.json               SOTA strategy comparison data   Table IV
saturate_fill_results.json          Weight-class saturation fill    SM Sec. III
floquet_dtc_results_ibm_kingston_20260401_113933.json
                                    DTC hardware XY correlators     SM Sec. XI
floquet_sim_sweep.json              DTC simulation sweep            SM Sec. XI
oq_n3_w_results.json                Rigetti n=3 W state (OQ)       Table I

============================================================
FORMAL VERIFICATION (verification/)
============================================================

proofs_registry.yaml       Central registry of all 11 formal manuscript claims,
                           keyed to manuscript line numbers and exact values.
                           Single source of truth for the six-tier pipeline.

certified_values.json      SageMath/FLINT exact arithmetic results (44 tests).
                           GMP integer and Pari number-field computations of
                           eigenvalues, N_4(K) counts, Gram matrix entries.

verification_results.json  Full six-tier pipeline verdict (~349 checks, zero
                           failures). Includes per-tier pass/fail, elapsed time,
                           platform metadata, and registry SHA-256.

lean4/                     Complete Lean 4 (Mathlib v4.29.0-rc8) source for all
                           11 formal statements (0 sorry). 17 files, 2236 lines.
                           Build: install elan, then `lake build` in lean4/.
                           Axioms.lean: 6 standard physical axioms (cited refs).
                           Statements.lean: consolidated formal claim list.

============================================================
AUXILIARY FILES (auxiliary/)
============================================================

Intermediate/exploratory results retained for reproducibility but not
part of the final submission. Includes: simulation outputs, pre-enrichment
raw data files, bootstrap analyses, checkpoint files, SOTA comparison data.

Pre-enrichment raw files (also present in top-level data/):
  sweep_ibm_fez_20260309_194458.json       Raw sweep before enrichment
  three_arm_hardware_ibm_fez_20260311_142656.json  Raw three-arm before enrichment
  w_repeat_results.json                    Raw W-repeat before enrichment
The enriched variants (with provenance blocks) are the authoritative versions.

============================================================
HARDWARE DATA SCHEMA
============================================================

Each single-run JSON file contains:
  backend              IBM backend name
  n_qubits             Number of qubits
  n_shots              Shots per circuit (1000)
  f_structured         Fidelity using structured measurement selection
  f_rand               Fidelity using random measurement selection
  advantage            f_structured - f_rand
  structured_expectations   Raw Pauli expectation values (structured set)
  rand_expectations         Raw Pauli expectation values (random set)
  timestamp            UTC timestamp of experiment
  provenance           Chain of custody (job IDs, raw counts, calibration)

The Bell pair file (220545) additionally contains:
  mode                 "full_tomography" (all 16 two-qubit Paulis measured)
  f_full               Fidelity from full tomography (all 16 Paulis)
  all_expectations     Expectation values for all 16 Pauli operators

The three-arm enriched file (20260311) contains per-seed results for
structured (K*), allocation-respecting random (AR), and uniform random
(UR) arms:
  runs[]               Per-seed fidelities for each arm (seeds 100, 200, 300)
  summary              Mean/std across seeds, allocation fraction

The Rigetti file (oq_grouped_results_20260316) contains same-session
basis-grouped results from Ankaa-3 via Open Quantum:
  results.n4_kstar     K* arm: fidelity, n_operators (137), n_bases (29/29)
  results.n4_rand      Random arm: fidelity, n_operators (137), n_bases (50/50)
  n4_kstar_raw_counts  Per-basis bitstring counts (29 bases x 1000 shots)
  n4_rand_raw_counts   Per-basis bitstring counts (50 bases x 1000 shots)
  compression          Basis compression ratios for both arms
  random_basis_selection  Seed=42 for random operator selection

The repeat-run enriched file (w_repeat_results_enriched) contains four
independent W-state runs with different random seeds (100, 200, 300, 400):
  runs[]               Per-run: seed, f_kstar, f_rand, delta_f, expectations
  summary              Population statistics (mean, std)
  provenance.raw_counts  Per-job bitstring histograms (12 jobs, 1000 shots each)
  random_basis_selection  Per-run seeds for random operator selection

============================================================
RANDOM SEED DOCUMENTATION
============================================================

All random operator selections are reproducible:

  W-repeat (IBM):     Per-run seeds [100, 200, 300, 400]
                      select_random_paulis(4, 137, seed=<run_seed>)
  Three-arm (IBM):    Per-run seeds [100, 200, 300]
                      select_random_paulis(4, 137, seed=<run_seed>)
  Rigetti 4q:         seed=42
                      select_random_paulis(4, 137, seed=42)
  Rigetti 8q:         seed=2027 (shuffle_seed=2026, random_seed=2027)
                      _generate_random_bases(115, 8, seed=2027)
  Sweep (IBM):        seed=42
  March 7 (IBM):      Exact seed not recorded (raw counts preserved)

============================================================
SOFTWARE VERSIONS
============================================================

Python 3.13.5, Qiskit 2.3.0, qiskit-ibm-runtime 0.45.1
NumPy 2.2.6, SciPy 1.15.2, SymPy 1.14.0, cvxpy 1.8.1
See scripts/requirements.txt for pinned dependencies.

============================================================
NOTES
============================================================

1. Fidelity values in hardware files are from maximum-likelihood
   reconstruction of raw expectation values. Table I in the paper
   reports fidelities from noise-robust MLE with shot-noise variance
   weighting, which may differ slightly from the raw values here.

2. The structured and random measurement sets each contain M Pauli
   operators (M=137 for 4-qubit states, M=9 for 2-qubit Bell).
   The complete 137-operator K* set is provided in kstar_operators_n4.json.

3. Both methods use identical circuit counts, shot counts, and
   reconstruction algorithms. The only variable is which Pauli
   operators are measured.

4. The compositional file (223638) reports per-patch fidelities from
   the 8-qubit scalability test (Sec. V C). One patch was not
   measured due to IBM runtime quota limits.

5. No readout error mitigation was applied to any hardware experiment.
   Both K* and random arms use identical raw measurement counts.

6. IBM 4q physical qubits: [0,1,2,3] (re-transpiled, optimization_level=1).
   IBM 8q physical qubits: [140,141,142,143,144,145,146,147].
   Rigetti 4q: [0,1,2,3]. Rigetti 8q: [0,1,2,3,4,5,6,7].

7. Rigetti/OQ does not expose backend calibration (T1, T2, gate errors)
   via API. This is a documented platform limitation.

8. IBM job results may expire after the retention period. All raw
   bitstring counts are preserved locally in the provenance blocks
   of each data file and will be archived on Zenodo.
