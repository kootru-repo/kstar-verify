/-
  KstarFormal.Combinatorics.CwClosedForm — Phase F: universal-n closed forms
  K* Verification: Krawtchouk spectral correspondence

  Phase F (post-E airtightness, 2026-04-17): symbolic closed forms
  for c_w_gen at K = 5, valid at all n:
    c_0(n, 5) = 1 + 2n
    c_1(n, 5) = 4n^2 - 2n  (written as 2n(2n - 1))
    c_2(n, 5) = 2n(n - 1)

  Each induction step walks the c_w_gen recursion down to the
  K-frontier {K - j^2 : j with j^2 <= K, j >= 1}, which for K = 5
  reaches K = 1 (via j = 2) and K = 4 (via j = 1).  We therefore
  need auxiliary closed forms at K in {0, 1, 3, 4} as well, which
  themselves reach K in {0, 3} (via j = 1) and K = 0 (via j = 2
  at K = 4).

  Strategy: for each (K, w) pair we need, state the recursion's
  unfold as a "succ-step" lemma whose RHS eliminates the Finset
  sums.  Since the filters at small K reduce to concrete singletons
  (or empty), a `decide` proves the filter equality, after which
  Finset.sum_singleton collapses the sum.

  All proofs: 0 sorry, 0 new propositional axioms.  The only axioms
  used are the Lean core + native_decide certificates from the
  decidable filter evaluations.
-/
import KstarFormal.Combinatorics.LatticeCountGeneric
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-! ## Filter-evaluation lemmas at each K needed by the recursion -/

/-- Empty even-j filter at K = 0 (range is {0}, fails j ≥ 1). -/
private theorem K0_even_filter :
    (Finset.range (Nat.sqrt 0 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 0 ∧ j % 2 = 0) = (∅ : Finset ℕ) := by
  native_decide

/-- Empty odd-j filter at K = 0. -/
private theorem K0_odd_filter :
    (Finset.range (Nat.sqrt 0 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 0 ∧ j % 2 = 1) = (∅ : Finset ℕ) := by
  native_decide

/-- Empty even-j filter at K = 1. -/
private theorem K1_even_filter :
    (Finset.range (Nat.sqrt 1 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 1 ∧ j % 2 = 0) = (∅ : Finset ℕ) := by
  native_decide

/-- Odd-j filter at K = 1 reduces to {1}. -/
private theorem K1_odd_filter :
    (Finset.range (Nat.sqrt 1 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 1 ∧ j % 2 = 1) = {1} := by
  native_decide

/-- Empty even-j filter at K = 3 (j = 2 has j*j = 4 > 3). -/
private theorem K3_even_filter :
    (Finset.range (Nat.sqrt 3 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 3 ∧ j % 2 = 0) = (∅ : Finset ℕ) := by
  native_decide

/-- Odd-j filter at K = 3 reduces to {1}. -/
private theorem K3_odd_filter :
    (Finset.range (Nat.sqrt 3 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 3 ∧ j % 2 = 1) = {1} := by
  native_decide

/-- Even-j filter at K = 4 reduces to {2}. -/
private theorem K4_even_filter :
    (Finset.range (Nat.sqrt 4 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 4 ∧ j % 2 = 0) = {2} := by
  native_decide

/-- Odd-j filter at K = 4 reduces to {1}. -/
private theorem K4_odd_filter :
    (Finset.range (Nat.sqrt 4 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 4 ∧ j % 2 = 1) = {1} := by
  native_decide

/-- Even-j filter at K = 5 reduces to {2}. -/
private theorem K5_even_filter :
    (Finset.range (Nat.sqrt 5 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 5 ∧ j % 2 = 0) = {2} := by
  native_decide

/-- Odd-j filter at K = 5 reduces to {1}. -/
private theorem K5_odd_filter :
    (Finset.range (Nat.sqrt 5 + 1)).filter
        (fun j => 1 ≤ j ∧ j * j ≤ 5 ∧ j % 2 = 1) = {1} := by
  native_decide

/-! ## Per-(K, w) recursion-step unfoldings

  Each lemma instantiates the c_w_gen recursion at a specific K and w
  shape, eliminating the Finset sum via the above filter lemmas.
-/

private theorem c_w_gen_succ_K0_w0 (n : ℕ) :
    c_w_gen (n + 1) 0 0 = c_w_gen n 0 0 := by
  simp only [c_w_gen, K0_even_filter, Finset.sum_empty, mul_zero, add_zero]

private theorem c_w_gen_succ_K0_wpos (n wp : ℕ) :
    c_w_gen (n + 1) 0 (wp + 1) = c_w_gen n 0 (wp + 1) := by
  simp only [c_w_gen, K0_even_filter, K0_odd_filter, Finset.sum_empty,
             mul_zero, add_zero]

private theorem c_w_gen_succ_K1_w0 (n : ℕ) :
    c_w_gen (n + 1) 1 0 = c_w_gen n 1 0 := by
  simp only [c_w_gen, K1_even_filter, Finset.sum_empty, mul_zero, add_zero]

private theorem c_w_gen_succ_K1_w1 (n : ℕ) :
    c_w_gen (n + 1) 1 1 = c_w_gen n 1 1 + 2 * c_w_gen n 0 0 := by
  simp only [c_w_gen, K1_even_filter, K1_odd_filter, Finset.sum_empty,
             Finset.sum_singleton, mul_zero, add_zero]

private theorem c_w_gen_succ_K1_wge2 (n w : ℕ) :
    c_w_gen (n + 1) 1 (w + 2) = c_w_gen n 1 (w + 2) + 2 * c_w_gen n 0 (w + 1) := by
  simp only [c_w_gen, K1_even_filter, K1_odd_filter, Finset.sum_empty,
             Finset.sum_singleton, mul_zero, add_zero]

private theorem c_w_gen_succ_K3_w0 (n : ℕ) :
    c_w_gen (n + 1) 3 0 = c_w_gen n 3 0 := by
  simp only [c_w_gen, K3_even_filter, Finset.sum_empty, mul_zero, add_zero]

private theorem c_w_gen_succ_K4_w0 (n : ℕ) :
    c_w_gen (n + 1) 4 0 = c_w_gen n 4 0 + 2 * c_w_gen n 0 0 := by
  simp only [c_w_gen, K4_even_filter, Finset.sum_singleton, mul_zero, add_zero]

private theorem c_w_gen_succ_K4_w1 (n : ℕ) :
    c_w_gen (n + 1) 4 1 = c_w_gen n 4 1 + 2 * c_w_gen n 0 1 + 2 * c_w_gen n 3 0 := by
  -- K4_w1 is the only step where every `+ 2*c_w_gen n _ _` summand is
  -- non-trivially present, so `mul_zero`/`add_zero` are genuinely unused
  -- here (the linter would flag them). Do NOT add them back for cosmetic
  -- consistency with the other `_succ_Kx_y` lemmas.
  simp only [c_w_gen, K4_even_filter, K4_odd_filter, Finset.sum_singleton]

private theorem c_w_gen_succ_K5_w0 (n : ℕ) :
    c_w_gen (n + 1) 5 0 = c_w_gen n 5 0 + 2 * c_w_gen n 1 0 := by
  simp only [c_w_gen, K5_even_filter, Finset.sum_singleton, mul_zero, add_zero]

private theorem c_w_gen_succ_K5_wpos (n wp : ℕ) :
    c_w_gen (n + 1) 5 (wp + 1) =
      c_w_gen n 5 (wp + 1) + 2 * c_w_gen n 1 (wp + 1)
        + 2 * c_w_gen n 4 wp := by
  simp only [c_w_gen, K5_even_filter, K5_odd_filter, Finset.sum_singleton]

/-! ## Closed-form theorems

  Now that the recursion-step unfoldings are explicit, each closed-form
  theorem reduces to a straightforward induction on n.
-/

theorem c_w_gen_n_K0_w0 (n : ℕ) : c_w_gen n 0 0 = 1 := by
  induction n with
  | zero => rfl
  | succ k ih => rw [c_w_gen_succ_K0_w0, ih]

theorem c_w_gen_n_K0_wpos (n w : ℕ) : c_w_gen n 0 (w + 1) = 0 := by
  induction n with
  | zero => rfl
  | succ k ih => rw [c_w_gen_succ_K0_wpos, ih]

theorem c_w_gen_n_K1_w0 (n : ℕ) : c_w_gen n 1 0 = 1 := by
  induction n with
  | zero => rfl
  | succ k ih => rw [c_w_gen_succ_K1_w0, ih]

theorem c_w_gen_n_K1_w1 (n : ℕ) : c_w_gen n 1 1 = 2 * n := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K1_w1, ih, c_w_gen_n_K0_w0]
      ring

theorem c_w_gen_n_K1_wge2 (n w : ℕ) : c_w_gen n 1 (w + 2) = 0 := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K1_wge2, ih, c_w_gen_n_K0_wpos]

theorem c_w_gen_n_K3_w0 (n : ℕ) : c_w_gen n 3 0 = 1 := by
  induction n with
  | zero => rfl
  | succ k ih => rw [c_w_gen_succ_K3_w0, ih]

theorem c_w_gen_n_K4_w0 (n : ℕ) : c_w_gen n 4 0 = 1 + 2 * n := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K4_w0, ih, c_w_gen_n_K0_w0]
      ring

theorem c_w_gen_n_K4_w1 (n : ℕ) : c_w_gen n 4 1 = 2 * n := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K4_w1, ih, c_w_gen_n_K0_wpos, c_w_gen_n_K3_w0]
      ring

/-! ## Phase F targets at K = 5 -/

/-- **Phase F target 1**: c_w_gen n 5 0 = 1 + 2n. -/
theorem c_w_gen_n_K5_w0 (n : ℕ) : c_w_gen n 5 0 = 1 + 2 * n := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K5_w0, ih, c_w_gen_n_K1_w0]
      ring

/-- **Phase F target 2**: c_w_gen n 5 1 = 2n(2n - 1). -/
theorem c_w_gen_n_K5_w1 (n : ℕ) : c_w_gen n 5 1 = 2 * n * (2 * n - 1) := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K5_wpos, ih, c_w_gen_n_K1_w1, c_w_gen_n_K4_w0]
      cases k with
      | zero => rfl
      | succ m =>
          have h1 : 2 * (m + 1) - 1 = 2 * m + 1 := by omega
          have h2 : 2 * (m + 1 + 1) - 1 = 2 * m + 3 := by omega
          rw [h1, h2]
          ring

/-- **Phase F target 3**: c_w_gen n 5 2 = 2n(n - 1). -/
theorem c_w_gen_n_K5_w2 (n : ℕ) : c_w_gen n 5 2 = 2 * n * (n - 1) := by
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [c_w_gen_succ_K5_wpos, ih, c_w_gen_n_K1_wge2, c_w_gen_n_K4_w1]
      cases k with
      | zero => rfl
      | succ m =>
          have h1 : m + 1 - 1 = m := by omega
          have h2 : m + 1 + 1 - 1 = m + 1 := by omega
          rw [h1, h2]
          ring

/-! ## Concrete n = 4 compatibility sanity check -/

theorem c_w_K5_closed_matches_n4 :
    c_w_gen 4 5 0 = 9 ∧ c_w_gen 4 5 1 = 56 ∧ c_w_gen 4 5 2 = 24 := by
  refine ⟨?_, ?_, ?_⟩
  · rw [c_w_gen_n_K5_w0]
  · rw [c_w_gen_n_K5_w1]
  · rw [c_w_gen_n_K5_w2]

/-! ## Phase F summary

  The three theorems `c_w_gen_n_K5_w0`, `c_w_gen_n_K5_w1`, and
  `c_w_gen_n_K5_w2` are the load-bearing Phase F deliverables:
  they replace the native_decide witnesses at n ∈ {4,5,6,7,8}
  in `UniversalKstar.lean` with universal-n closed forms.

  The greedy-redistribution saturation for w ∈ {1, 2} at K = 5 is
  the arithmetic consequence of these closed forms (after the
  polynomial inequality 3n^2 - n ≥ 0 at n ≥ 1). That lift to
  explicit greedy output happens in `UniversalKstar.lean` where
  `greedyRedist` is computed over concrete lists; Phase F supplies
  the underlying closed forms, which are universal over n.
-/
