# Verification Coverage Matrix

Maps every formal statement and empirical claim in the manuscript to its
machine-verification tests. No claim is verified by AI — all checks are
deterministic code runnable by any reviewer.

## Formal statements

| Statement | Proof scope | Tier 0 (Lean4) | Tier 1 (SageMath) | Tier 2 (SymPy) | Tier 3 (NumPy) | Tier 4 (Indep.) | Tier 5 (HS→F) | Tier 6 (Chain) |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| **Lemma 1** (Fisher-Hessian) | general n | pauli_trace_offdiag | hessian_diagonal | — | test_hessian (11) | verify_fisher_frame | — | — |
| **Theorem 1(i)** (basin sep.) | axiom-based | thm1_i_informative_identifiability | condition_number, W_high_weight | test_basin_separation (21) | test_hessian (flat dirs) | verify_operator_set | — | chain2 |
| **Theorem 1(ii)** (random miss) | concrete n=4 | thm1_ii_expected_missing | hypergeometric_exact, hypergeometric_bound | test_basin_separation | — | — | — | chain3 |
| **Theorem 1(iii)** (HS error) | axiom-based | thm1_iii_eps_pos_chain | purity_upper, approx_local | test_purity_bound (12) | test_fidelity_bounds (14) | — | **all 10 checks** | chain2 |
| **Corollary 1** (approx. local) | axiom-based | cor1_approx_locality | approx_local_bound, eps_tail | test_purity_bound | test_fidelity_bounds | — | HS-Pauli identity | chain2 step5 |
| **Lemma 2(i)** (purity bound) | general n | purity_bound_prob | purity_upper, depol_eps_pos | test_purity_bound | test_fidelity_bounds | — | near-pure recovery | — |
| **Lemma 2(ii)** (tightness) | concrete n=4 | — | pure_state_completeness | test_purity_bound | test_fidelity_bounds (adversarial) | — | — | — |
| **Lemma 3** (monotonicity) | concrete n=4 | eigenvalue_monotone_K4_K5 | eigenvalue_monotone, greedy_monotone | — | test_linear_bound (monotone) | — | — | — |
| **Theorem 2(i)** (mu_w<1 → F≤½) | axiom-based | witness_fidelity_zero | — | — | test_linear_bound (Pauli-axis) | — | — | — |
| **Theorem 2(ii)** (mu_w=1 → HS bound) | axiom-based | no_estimator_above_half | spectral_hs_bound | — | test_linear_bound | — | — | — |
| **Theorem 2(iii)** (phase transition) | concrete n=4 | spectral_sharp_gap | — | — | test_linear_bound (sharp gap) | — | — | — |
| **Theorem 2(iv)** (K* phases) | concrete n=4 | weight_saturation_K5 | phase2_values, phase3_saturation | test_weight_saturation (10) | test_linear_bound (Delsarte) | — | — | chain1 |
| **Lemma 4** (coupon-collector) | general n | hyper_product_bound | hypergeometric_exact, random_vanishing | test_asymptotic | — | — | — | chain3 |
| **Corollary 2** (lower bound) | concrete n=4 | operator_lower_bound_66 | operator_lower_bound | test_asymptotic | — | — | — | — |
| **Theorem 3** (asymptotic) | concrete n=4 | saturation_ratio_lt_one | mn_over_n_ratio | test_asymptotic | — | — | — | chain3 |
| **Lemma 5** (general q) | concrete n=4 | — | krawtchouk_orthogonality, spectral_mass | test_universality (21) | test_linear_bound (q=3,5) | verify_qutrit | — | — |
| **Lemma 6** (support-complete) | concrete n=4 | gram_rank_16_K5 | kstar_full_rank_cross | test_gram_krawtchouk (12) | — | — | — | — |

### Proof scope legend

- **concrete n=4**: Lean4 uses `native_decide` to verify the claim at n=4 (the paper's operational dimension). General-n is verified numerically by Tiers 2-3.
- **general n**: Lean4 proves for arbitrary n (structural proof, not computation).
- **axiom-based**: Lean4 proves algebraic premises; quantum conclusion uses declared axioms (Uhlmann, Fuchs-van de Graaf, etc.). See `lean4/KstarFormal/Axioms.lean`.

## Key numerical values

| Value | Manuscript | Tier 1 | Tier 2 | Tier 4 |
|-------|-----------|:-:|:-:|:-:|
| N₄(5) = 137 | Table I | N4_K5 | check 13 | verify_operator_set |
| M_w = [1,12,54,54,16] | Eq. 12 | phase3_M_w | check 67 | verify_operator_set |
| c_w = [9,56,24,32,16] | Appendix A | phase3_saturation | check 11, 73 | verify_operator_set |
| A_w = [1,12,54,108,81] | Sec. II | operator_lower_bound | check 14, 68 | verify_operator_set |
| λ_w = [144,224,64,128,256] | L803 (formula) | eigenvalues_n4_K5 | check 9 | verify_propositions |
| κ_G = 4 | Sec. III B | gram_condition_number | — | — |
| w_sat = 2 | Theorem 1 | kstar_minimality | check 21, 72 | verify_operator_set |
| 29 K* bases | Sec. II | — | — | verify_operator_set |
| ε_pos(W,k=2) = 11/256 | Lemma 2 | purity_upper | check 39-40 | — |
| ε_pos(product,k=2) = 5/256 | Lemma 2 | purity_upper | check 37-38 | — |
| ε_pos(GHZ,k=2) = 9/256 | Lemma 2 | purity_upper | check 41-42 | — |
| P(all wt-1) < 5×10⁻⁴ | Theorem 1(ii) | hypergeometric_exact | — | verify_fisher_frame |
| 56 informative W ops | Lemma 1 | condition_number | — | — |
| κ_info(W) = 1 | Lemma 1 | condition_number | — | — |
| K*=q² (q≥3) | Lemma 5 | spectral_mass_threshold | check 46-49 | verify_qutrit |
| 2d=2^{d-1} → d=4 | Thm rigidity | dimensional_rigidity | — | verify_propositions |

## Hardware / empirical claims

| Claim | Manuscript | Tier 4 (Indep.) | Tier 6 (Chain) |
|-------|-----------|:-:|:-:|
| F(K*,W) = 0.872±0.021 | Table III | verify_hardware_fidelities | chain4 (W-repeat) |
| F(rand,W) = 0.542±0.076 | Table III | verify_hardware_fidelities | chain4 (W-repeat) |
| ΔF(W) = +0.33±0.07 | Table III | verify_hardware_fidelities | chain4 (dF check) |
| F(K*,product) > 0.99 | Table III | verify_hardware_fidelities | — |
| 8q compositional F > 0.99 | Sec. IV C | verify_hardware_fidelities | — |
| F(K*) > F(AR) > F(UR) | Three-arm | verify_allocation | chain4 (three-arm) |
| Allocation fraction ~76% | Three-arm | verify_allocation | chain4 |
| GHZ F ≈ 0.50 (anomaly) | Sec. IV B | verify_ghz_reanalysis | — |
| RDM fidelities (1-/2-RDM) | SM Sec. III | verify_rdm_fidelities | — |
| SOTA 7-strategy ranking | Table II | verify_sota_and_stats | — |
| Rigetti F(K*)=0.816 | Sec. IV | verify_hardware_fidelities | — |
| Figure data consistency | All figures | verify_figures | — |

## Supplemental material

| SM Section | Content | Verified by |
|-----------|---------|-------------|
| S1: Per-weight budget | Allocation derivation | Tier 2: test_weight_saturation; Tier 4: verify_operator_set |
| S2: Support-completeness | Lemma 6 | Tier 1: kstar_full_rank_cross; Tier 2: test_gram_krawtchouk |
| S3: Three layers of K* | K*_full vs K*_paper | Tier 1: kstar_minimality; Tier 2: test_weight_saturation |
| S4: Circuit compression | 137→29 bases | Tier 4: verify_operator_set (greedy_basis_cover) |
| S5: Krawtchouk tables | Polynomial values | Tier 1: krawtchouk tests; Tier 2: test_krawtchouk |
| S6: RDM fidelities | 1-/2-RDM | Tier 4: verify_rdm_fidelities |
| S7: Fisher/frame theory | Frame potential, coherence | Tier 4: verify_fisher_frame |
| Extended tables | Ququint, scaling, ordering | Tier 4: verify_sm_extended |

## Residual gaps (not machine-verifiable)

| Gap | Why | Mitigation |
|-----|-----|------------|
| Proof text correctness | Natural-language arguments require human reading | proofs.tex + remediation.md document all 13 audited issues |
| QMA-completeness of PSD feasibility | Lemma 2(ii) tightness direction — no polynomial algorithm exists | Acknowledged in theorem statement; upper bound (the useful direction) is unconditional |
| M_n/N ≤ c < 1 for all n | Theorem 3(ii) hypothesis — numerically verified n=3..6, not proven | Explicitly stated as hypothesis; verified in Tier 6 chain3 for n=3..6 |
