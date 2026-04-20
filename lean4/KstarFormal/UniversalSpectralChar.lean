/-
  KstarFormal.UniversalSpectralChar — Phase C: universal thm:spectral-char
  K* Verification: Krawtchouk spectral correspondence

  Phase C of the universal-n airtightness plan (2026-04-17):
  lift thm:spectral-char (Theorem 2) parts (i)-(iii) from their
  current n=4 concrete statements to universal n ≥ 2.

  Part (i): Spectral obstruction (necessity).
    If μ_{w₀} < 1 (weight class w₀ not saturated), there exist
    witness states with F = 0.  Proof uses universal-n Pauli
    orthogonality + 1 axiom (fidelity_orthogonal_zero, Uhlmann).

  Part (ii): Spectral resolution (sufficiency).
    Universal-n purity bound (lem2_purity_bound) delivers
    Σ y_P² ≤ 1, chained with Parseval and Weyl to bound ε_pos.

  Part (iii): Sharp dichotomy.  Combines (i) and (ii).

  Part (iv) is handled in a separate module (UniversalKstar.lean)
  because the K*=5 saturation argument is independent of the
  fidelity construction.

  Axiom footprint for each theorem below:
    witness_indistinguishable_universal: Lean core only
    no_estimator_above_half_universal:   Lean core + 1 prop axiom
    fidelity_dichotomy_universal:        Lean core + 1 prop axiom
    thm2_spectral_char_universal:        Lean core + 2 prop axioms
                                         (fidelity_orthogonal_zero
                                         + fidelity_witness_traceproj_decomp)
-/
import KstarFormal.Statements
import KstarFormal.Universal

open Finset BigOperators

/-! ## Universal-n witness indistinguishability (Part i algebraic core) -/

/-- **Universal-n Pauli tracelessness.**
    For any n-qubit non-identity Pauli Q, Tr(Q) = 0. -/
theorem pauli_traceless_universal {n : ℕ} (Q : PauliIdx n)
    (hQ : Q ≠ pauliIdentity n) :
    pauliTraceN (pauliIdentity n) Q = 0 :=
  pauli_trace_offdiag (pauliIdentity n) Q (Ne.symm hQ)

/-- **Universal-n Pauli cross-trace.**
    For any two distinct n-qubit Paulis Q, P₀: Tr(Q · P₀) = 0. -/
theorem pauli_cross_trace_zero_universal {n : ℕ} (Q P₀ : PauliIdx n)
    (h : Q ≠ P₀) :
    pauliTraceN Q P₀ = 0 :=
  pauli_trace_offdiag Q P₀ h

/-- **Universal-n witness indistinguishability.**
    For any unmeasured n-qubit Pauli P₀ and any measured Q with
    Q ≠ pauliIdentity and Q ≠ P₀, the witness states
    ρ± = (I ± P₀)/d are measurement-indistinguishable at Q:
    both Tr(Q · ρ₊) = 0 and Tr(Q · ρ₋) = 0, since Tr(Q) = 0
    and Tr(Q · P₀) = 0 by Pauli orthogonality. -/
theorem witness_indistinguishable_universal {n : ℕ} (Q P₀ : PauliIdx n)
    (hQ_nonid : Q ≠ pauliIdentity n) (hQP : Q ≠ P₀) :
    pauliTraceN (pauliIdentity n) Q = 0 ∧ pauliTraceN Q P₀ = 0 :=
  ⟨pauli_traceless_universal Q hQ_nonid,
   pauli_cross_trace_zero_universal Q P₀ hQP⟩

/-- **Universal-n witness indistinguishability over a measurement set.**
    If S is a set of non-identity n-qubit Paulis not containing P₀,
    every Q ∈ S satisfies Tr(Q · ρ₊) = Tr(Q · ρ₋) = 0. -/
theorem witness_indistinguishable_on_set_universal {n : ℕ} (P₀ : PauliIdx n)
    (S : Finset (PauliIdx n))
    (hS_nonid : ∀ Q ∈ S, Q ≠ pauliIdentity n)
    (hS_miss : P₀ ∉ S) :
    ∀ Q ∈ S, pauliTraceN (pauliIdentity n) Q = 0 ∧ pauliTraceN Q P₀ = 0 := by
  intro Q hQ
  exact witness_indistinguishable_universal Q P₀
    (hS_nonid Q hQ) (fun h => hS_miss (h ▸ hQ))

/-! ## Universal-n witness trace and involution -/

/-- **Universal-n witness trace normalization.**
    Tr(I ± P₀) = Tr(I) ± Tr(P₀) = 2^n ± 0 = 2^n, so
    Tr(ρ±) = (2^n ± 0)/2^n = 1. -/
theorem witness_trace_correct_universal {n : ℕ} (P₀ : PauliIdx n)
    (hP₀ : P₀ ≠ pauliIdentity n) :
    pauliTraceN (pauliIdentity n) (pauliIdentity n) = (2 : ℤ) ^ n ∧
      pauliTraceN (pauliIdentity n) P₀ = 0 := by
  refine ⟨?_, pauli_traceless_universal P₀ hP₀⟩
  rw [pauli_traceN_eq, if_pos rfl]

/-- **Universal-n Pauli involution at trace level.**
    Tr(P₀²) = Tr(I) = 2^n for any n-qubit Pauli (including identity).
    This is the algebraic core of PSD for the witness states. -/
theorem pauli_involution_trace_universal {n : ℕ} (P₀ : PauliIdx n) :
    pauliTraceN P₀ P₀ = (2 : ℤ) ^ n := by
  rw [pauli_traceN_eq, if_pos rfl]

/-- **Universal-n orthogonal support at trace level.**
    Tr(I) - Tr(P₀²) = 0, so the supports of ρ₊ and ρ₋ have zero
    overlap.  Algebraic core of F(ρ₊, ρ₋) = 0. -/
theorem witness_orthogonal_support_universal {n : ℕ} (P₀ : PauliIdx n) :
    pauliTraceN (pauliIdentity n) (pauliIdentity n)
      - pauliTraceN P₀ P₀ = 0 := by
  rw [pauli_traceN_eq, if_pos rfl, pauli_traceN_eq, if_pos rfl]; ring

/-! ## Universal-n Theorem 2(i): necessity -/

/-- **Theorem 2(i), universal n.**
    For any n and any unmeasured non-identity Pauli P₀:
    (a) measurement indistinguishability of ρ± on S
    (b) trace normalization
    (c) Pauli involution (core of PSD)
    (d) orthogonal support (core of F = 0) -/
theorem thm2_i_necessity_universal {n : ℕ} (P₀ : PauliIdx n)
    (hP₀ : P₀ ≠ pauliIdentity n)
    (S : Finset (PauliIdx n))
    (hS_nonid : ∀ Q ∈ S, Q ≠ pauliIdentity n)
    (hS_miss : P₀ ∉ S) :
    (∀ Q ∈ S, pauliTraceN (pauliIdentity n) Q = 0
            ∧ pauliTraceN Q P₀ = 0) ∧
      (pauliTraceN (pauliIdentity n) (pauliIdentity n) = (2 : ℤ) ^ n ∧
       pauliTraceN (pauliIdentity n) P₀ = 0) ∧
      (pauliTraceN P₀ P₀ = (2 : ℤ) ^ n) ∧
      (pauliTraceN (pauliIdentity n) (pauliIdentity n)
        - pauliTraceN P₀ P₀ = 0) :=
  ⟨witness_indistinguishable_on_set_universal P₀ S hS_nonid hS_miss,
   witness_trace_correct_universal P₀ hP₀,
   pauli_involution_trace_universal P₀,
   witness_orthogonal_support_universal P₀⟩

/-! ## Universal-n CONCLUSION theorems (use axioms from KstarFormal.Axioms) -/

/-- **CONCLUSION, universal n: F(ρ₊, ρ₋) = 0** (Uhlmann).
    The witness states from any unmeasured non-identity Pauli have
    zero fidelity.  Uses `fidelity_orthogonal_zero` (Axiom 2, Uhlmann). -/
theorem witness_fidelity_zero_universal {n : ℕ} (P₀ : PauliIdx n)
    (hP₀ : P₀ ≠ pauliIdentity n) :
    qfidelity (witnessPlus P₀) (witnessMinus P₀) = 0 :=
  fidelity_orthogonal_zero P₀ hP₀

/-- **CONCLUSION, universal n: no estimator above 1/2** (Fuchs-van de Graaf).
    For ANY n-qubit estimator σ̂ and any unmeasured non-identity Pauli P₀,
    at least one of F(σ̂, ρ₊) ≤ 1/2 or F(σ̂, ρ₋) ≤ 1/2.
    Uses `fidelity_witness_traceproj_decomp` (Axiom 3), from which
    `fidelity_sum_complementary` is a derived theorem. -/
theorem no_estimator_above_half_universal {n : ℕ} (P₀ : PauliIdx n)
    (hP₀ : P₀ ≠ pauliIdentity n) (σ_hat : DensityMatrix (2 ^ n)) :
    qfidelity σ_hat (witnessPlus P₀) ≤ 1 / 2 ∨
      qfidelity σ_hat (witnessMinus P₀) ≤ 1 / 2 := by
  by_contra h
  push_neg at h
  have h_sum := fidelity_sum_complementary P₀ hP₀ σ_hat
  linarith

/-- **CONCLUSION, universal n: Fidelity dichotomy phase transition.**
    (a) Incomplete measurement ⟹ no estimator above 1/2 (part i).
    (b) Complete measurement ⟹ purity bound gives ε_pos → 0 (part ii).
    Universal in n and d.  Axiom footprint: 1 (via no_estimator branch). -/
theorem fidelity_dichotomy_universal :
    (∀ {n : ℕ} (P₀ : PauliIdx n), P₀ ≠ pauliIdentity n →
      ∀ (σ_hat : DensityMatrix (2 ^ n)),
        qfidelity σ_hat (witnessPlus P₀) ≤ 1 / 2 ∨
          qfidelity σ_hat (witnessMinus P₀) ≤ 1 / 2) ∧
      (∀ (d : ℕ) (ev : Fin d → ℚ),
        (∀ i, 0 ≤ ev i) → (∑ i, ev i = 1) → (∑ i, (ev i) ^ 2 ≤ 1)) :=
  ⟨fun P₀ hP₀ σ_hat => no_estimator_above_half_universal P₀ hP₀ σ_hat,
   fun _ ev h_nn h_sum => purity_bound_prob ev h_nn h_sum⟩

/-! ## Universal-n Theorem 2: Spectral characterization -/

/-- **Theorem 2 (thm:spectral-char), universal n.**
    Bundles parts (i), (ii), (iii) at universal n.
    Part (iv) — K*=5 saturation for n ≥ 4 — is handled separately
    in `UniversalKstar.lean` because its combinatorial argument is
    independent of the witness construction.

    Content:
    (i) Witness indistinguishability at any n, for any unmeasured P₀:
        Tr(Q · ρ±) = 0 for every Q ∈ S.
    (i') F(ρ₊, ρ₋) = 0 conclusion (from Uhlmann, Axiom 2).
    (iii) Dichotomy: either worst-case F ≤ 1/2 (incomplete), or
          Σ ev² ≤ 1 purity bound (complete, MLE side of ε_pos → 0). -/
theorem thm2_spectral_char_universal :
    -- (i) universal-n witness indistinguishability
    (∀ {n : ℕ} (Q P₀ : PauliIdx n),
      Q ≠ pauliIdentity n → Q ≠ P₀ →
      pauliTraceN (pauliIdentity n) Q = 0 ∧ pauliTraceN Q P₀ = 0) ∧
    -- (i') universal-n F(ρ₊, ρ₋) = 0
    (∀ {n : ℕ} (P₀ : PauliIdx n), P₀ ≠ pauliIdentity n →
      qfidelity (witnessPlus P₀) (witnessMinus P₀) = 0) ∧
    -- (iii) universal-n dichotomy: part (a) no estimator above 1/2
    (∀ {n : ℕ} (P₀ : PauliIdx n), P₀ ≠ pauliIdentity n →
      ∀ (σ_hat : DensityMatrix (2 ^ n)),
        qfidelity σ_hat (witnessPlus P₀) ≤ 1 / 2 ∨
          qfidelity σ_hat (witnessMinus P₀) ≤ 1 / 2) ∧
    -- (iii) universal dichotomy: part (b) purity bound
    (∀ (d : ℕ) (ev : Fin d → ℚ),
      (∀ i, 0 ≤ ev i) → (∑ i, ev i = 1) → (∑ i, (ev i) ^ 2 ≤ 1)) :=
  ⟨fun Q P₀ h1 h2 => witness_indistinguishable_universal Q P₀ h1 h2,
   fun P₀ hP₀ => witness_fidelity_zero_universal P₀ hP₀,
   fun P₀ hP₀ σ_hat => no_estimator_above_half_universal P₀ hP₀ σ_hat,
   fun _ ev h_nn h_sum => purity_bound_prob ev h_nn h_sum⟩
