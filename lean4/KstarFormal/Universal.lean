/-
  KstarFormal.Universal — Universal-n airtightness wrappers (Phase B)
  K* Verification: Krawtchouk spectral correspondence

  Phase B (universal-n airtightness plan, 2026-04-17) deliverables:
  explicit universal-n / universal-d / universal-K companions for the
  manuscript theorems whose paper-scope claim is universal.

  Most underlying proofs are already universal in the existing Lean
  sources; this file collects them under explicit `*_universal` names
  that match the paper's scope declarations, and adds the
  universal-K monotonicity proof for Lemma 3 (eigenvalue monotonicity).

  Each theorem's #print axioms output is expected to contain only
  { propext, Classical.choice, Quot.sound } plus (optionally) the
  existing propositional axioms in KstarFormal.Axioms. No new axioms.
-/
import KstarFormal.Statements
import KstarFormal.Combinatorics.LatticeCountGeneric
import Mathlib.Tactic.Linarith

open Finset

/-! ## Lemma 1 (lem:hessian) at universal n -/

/-- **Lemma 1 (lem:hessian), universal n.**
    For any n, the trace product of distinct n-qubit Paulis is zero:
    Tr(P · Q) = 0 whenever P ≠ Q. -/
theorem lem1_hessian_universal :
    ∀ (n : ℕ) (P Q : PauliIdx n), P ≠ Q → pauliTraceN P Q = 0 :=
  lem1_hessian_diagonal

/-! ## Lemma 2 (lem:purity_main) at universal d -/

/-- **Lemma 2 (lem:purity_main), universal d.** -/
theorem lem2_purity_universal :
    ∀ (d : ℕ) (ev : Fin d → ℚ),
      (∀ i, 0 ≤ ev i) → (∑ i, ev i = 1) → (∑ i, (ev i) ^ 2 ≤ 1) :=
  lem2_purity_bound

/-! ## Corollary 1 (cor:approx_local) at universal d -/

/-- **Corollary 1 (cor:approx_local), universal d.** -/
theorem cor1_approx_locality_universal :
    ∀ (d : ℕ), 0 < d →
    ∀ (tr_sq S_k unmeasured eps_pos : ℚ),
    tr_sq ≤ 1 →
    S_k + unmeasured = (d : ℚ) * tr_sq - 1 →
    eps_pos ≤ unmeasured →
    eps_pos ≤ (d : ℚ) - 1 - S_k :=
  cor1_approx_locality_stmt

/-! ## Theorem 1 (thm:basin) at universal n, d -/

/-- **Theorem 1 (thm:basin), universal parts.**
    Bundles parts (i), (ii), (iii) at their universal scopes. -/
theorem thm1_basin_universal :
    -- Part (i): universal Hessian diagonality
    (∀ (n : ℕ) (P Q : PauliIdx n), P ≠ Q → pauliTraceN P Q = 0) ∧
    -- Part (ii): universal expected-missing count
    (∀ (A_w M N : ℕ), M ≤ N → 0 < N →
      expected_missing_count A_w M N = A_w * ((N - M : ℤ) / N)) ∧
    -- Part (iii): universal HS finite-sample bound
    (∀ (d : ℕ), 0 < d →
      ∀ (tr_sq S_k unmeasured eps_pos : ℚ),
      tr_sq ≤ 1 →
      S_k + unmeasured = (d : ℚ) * tr_sq - 1 →
      eps_pos ≤ unmeasured →
      eps_pos ≤ (d : ℚ) - 1 - S_k) :=
  ⟨lem1_hessian_diagonal, thm1_ii_expected_missing_abstract, thm1_iii_eps_pos_chain⟩

/-! ## Lemma 3 (lem:monotone) at universal n, K

  The paper claims λ_w(K) is non-decreasing in K for every weight class w.
  This is lifted from the K=4→K=5 at n=4 Lean witness to a generic
  statement over all n, K, K'.

  Proof structure:
    1. c_w_gen monotonicity in K at fixed n, w (induction on n).
       The inductive step handles the recursion's three components:
         (A) the j = 0 term (direct IH on the first dimension),
         (B) the even-j shell sum (two-step calc: IH on argument,
             then subset inclusion of filter),
         (C) the odd-j shell sum (same pattern, with parity-shift).
    2. Lift to krawtchouk_eigenvalue_gen via the divisor-preserving
       monotonicity of Nat.div.
-/

/-- Helper: even-j shell filter on range (Nat.sqrt K + 1).
    `abbrev` so simp unfolds and re-folds transparently. -/
abbrev even_j_filter (K : ℕ) : Finset ℕ :=
  (Finset.range (Nat.sqrt K + 1)).filter (fun j => 1 ≤ j ∧ j * j ≤ K ∧ j % 2 = 0)

/-- Helper: odd-j shell filter on range (Nat.sqrt K + 1). -/
abbrev odd_j_filter (K : ℕ) : Finset ℕ :=
  (Finset.range (Nat.sqrt K + 1)).filter (fun j => 1 ≤ j ∧ j * j ≤ K ∧ j % 2 = 1)

/-- The even-j filter is monotone in K. -/
private theorem even_j_filter_subset {K K' : ℕ} (hK : K ≤ K') :
    even_j_filter K ⊆ even_j_filter K' := by
  intro j hj
  simp only [even_j_filter, Finset.mem_filter, Finset.mem_range] at hj ⊢
  refine ⟨?_, hj.2.1, le_trans hj.2.2.1 hK, hj.2.2.2⟩
  exact Nat.lt_succ_of_le (le_trans (Nat.le_of_lt_succ hj.1) (Nat.sqrt_le_sqrt hK))

/-- The odd-j filter is monotone in K. -/
private theorem odd_j_filter_subset {K K' : ℕ} (hK : K ≤ K') :
    odd_j_filter K ⊆ odd_j_filter K' := by
  intro j hj
  simp only [odd_j_filter, Finset.mem_filter, Finset.mem_range] at hj ⊢
  refine ⟨?_, hj.2.1, le_trans hj.2.2.1 hK, hj.2.2.2⟩
  exact Nat.lt_succ_of_le (le_trans (Nat.le_of_lt_succ hj.1) (Nat.sqrt_le_sqrt hK))

/-- Monotonicity of the generic parity-weight count in the ball radius K,
    for fixed dimension n and weight w. Proved by induction on n. -/
theorem c_w_gen_mono_K :
    ∀ (n K K' w : ℕ), K ≤ K' → c_w_gen n K w ≤ c_w_gen n K' w := by
  intro n
  induction n with
  | zero =>
      intro K K' w _
      cases w with
      | zero => simp [c_w_gen]
      | succ _ => simp [c_w_gen]
  | succ n ih =>
      intro K K' w hK
      have h_even :
          (even_j_filter K).sum (fun j => c_w_gen n (K - j * j) w)
            ≤ (even_j_filter K').sum (fun j => c_w_gen n (K' - j * j) w) := by
        calc (even_j_filter K).sum (fun j => c_w_gen n (K - j * j) w)
            ≤ (even_j_filter K).sum (fun j => c_w_gen n (K' - j * j) w) := by
              apply Finset.sum_le_sum
              intro j _
              exact ih _ _ _ (Nat.sub_le_sub_right hK _)
          _ ≤ (even_j_filter K').sum (fun j => c_w_gen n (K' - j * j) w) := by
              apply Finset.sum_le_sum_of_subset_of_nonneg (even_j_filter_subset hK)
              intros; exact Nat.zero_le _
      have h_odd : ∀ wp,
          (odd_j_filter K).sum (fun j => c_w_gen n (K - j * j) wp)
            ≤ (odd_j_filter K').sum (fun j => c_w_gen n (K' - j * j) wp) := by
        intro wp
        calc (odd_j_filter K).sum (fun j => c_w_gen n (K - j * j) wp)
            ≤ (odd_j_filter K).sum (fun j => c_w_gen n (K' - j * j) wp) := by
              apply Finset.sum_le_sum
              intro j _
              exact ih _ _ _ (Nat.sub_le_sub_right hK _)
          _ ≤ (odd_j_filter K').sum (fun j => c_w_gen n (K' - j * j) wp) := by
              apply Finset.sum_le_sum_of_subset_of_nonneg (odd_j_filter_subset hK)
              intros; exact Nat.zero_le _
      -- Unfold both sides' recursion and combine. Since even_j_filter /
      -- odd_j_filter are abbrevs, simp only [c_w_gen] reveals them in
      -- the unfolded form that h_even / h_odd match against.
      show c_w_gen (n+1) K w ≤ c_w_gen (n+1) K' w
      simp only [c_w_gen]
      refine Nat.add_le_add (Nat.add_le_add (ih K K' w hK) ?_) ?_
      · exact Nat.mul_le_mul_left 2 h_even
      · refine Nat.mul_le_mul_left 2 ?_
        cases w with
        | zero => simp
        | succ wp => exact h_odd wp

/-- **Lemma 3 (lem:monotone), universal.**
    The generic Krawtchouk eigenvalue is non-decreasing in K for any
    dimension n and any weight class w. Follows from c_w_gen_mono_K
    and the monotonicity of Nat.div. -/
theorem lem3_eigenvalue_monotonicity_universal :
    ∀ (n K K' w : ℕ), K ≤ K' →
      krawtchouk_eigenvalue_gen n K w ≤ krawtchouk_eigenvalue_gen n K' w := by
  intro n K K' w hK
  unfold krawtchouk_eigenvalue_gen
  exact Nat.div_le_div_right (Nat.mul_le_mul_left _ (c_w_gen_mono_K n K K' w hK))
