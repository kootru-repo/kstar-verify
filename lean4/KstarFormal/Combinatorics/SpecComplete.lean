/-
  KstarFormal.Combinatorics.SpecComplete — Spectral completeness threshold
  K* Verification: Krawtchouk spectral correspondence
  Registry: prop:spec-complete (SM Proposition 2)

  Formalizes the manuscript claim that for any q ≥ 2 and n ≥ 1, the
  Gram matrix G(K) of H(n,q) achieves spectral completeness
  (c_w(K) > 0 for all w = 0,...,n) at K = n, and that this threshold
  is tight: c_n(n-1) = 0.

  Both halves are elementary statements about integer lattice vectors:

    Sufficiency: for every 0 ≤ w ≤ n there exists m ∈ ℤⁿ with
                 |m|² ≤ n and parity weight w(m̄) = w.
                 Witness: m = (1,…,1,0,…,0) with w leading ones.

    Tightness:   any m ∈ ℤⁿ with parity weight n has all coordinates
                 nonzero mod q, hence each |mᵢ| ≥ 1, hence |m|² ≥ n.

  Status: Tier 1 formalization, sorry-free.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Data.Int.GCD
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Fintype.Card
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

namespace KstarFormal.SpecComplete

open Finset BigOperators

/-- Squared Euclidean norm of an integer vector indexed by `Fin n`. -/
def normSq {n : ℕ} (m : Fin n → ℤ) : ℤ := ∑ i, (m i) ^ 2

/-- Parity weight w_q(m): number of coordinates not divisible by q. -/
def parityWeight {n : ℕ} (q : ℕ) (m : Fin n → ℤ) : ℕ :=
  (Finset.univ.filter fun i => ¬ ((q : ℤ) ∣ m i)).card

/-- Canonical sufficiency witness: w leading ones, then n - w trailing zeros. -/
def witnessVec (n w : ℕ) : Fin n → ℤ :=
  fun i => if i.val < w then 1 else 0

/-! ## Counting `{i : Fin n // i.val < w}` -/

/-- For w ≤ n, the set of indices i ∈ Fin n with i.val < w has cardinality w. -/
lemma card_lt_eq (n w : ℕ) (hw : w ≤ n) :
    (Finset.univ.filter (fun i : Fin n => i.val < w)).card = w := by
  rw [← Fintype.card_subtype]
  let e : {i : Fin n // i.val < w} ≃ Fin w :=
    { toFun := fun ⟨i, h⟩ => ⟨i.val, h⟩
      invFun := fun ⟨k, h⟩ => ⟨⟨k, h.trans_le hw⟩, h⟩
      left_inv := fun _ => rfl
      right_inv := fun _ => rfl }
  rw [Fintype.card_congr e, Fintype.card_fin]

/-! ## Witness evaluation -/

/-- The squared norm of the canonical witness equals w (when w ≤ n). -/
lemma normSq_witness (n w : ℕ) (hw : w ≤ n) :
    normSq (witnessVec n w) = (w : ℤ) := by
  unfold normSq witnessVec
  have h1 : ∀ i : Fin n,
      (if i.val < w then (1 : ℤ) else 0) ^ 2
        = if i.val < w then (1 : ℤ) else 0 := by
    intro i
    by_cases h : i.val < w
    · simp [h]
    · simp [h]
  rw [Finset.sum_congr rfl (fun i _ => h1 i)]
  rw [← Finset.sum_filter]
  rw [Finset.sum_const]
  rw [card_lt_eq n w hw]
  simp

/-- The parity weight of the canonical witness equals w (when q ≥ 2 and w ≤ n). -/
lemma parityWeight_witness (q n w : ℕ) (hq : 2 ≤ q) (hw : w ≤ n) :
    parityWeight q (witnessVec n w) = w := by
  unfold parityWeight witnessVec
  have heq :
      (Finset.univ.filter fun i : Fin n =>
          ¬ ((q : ℤ) ∣ if i.val < w then (1 : ℤ) else 0))
        = Finset.univ.filter fun i : Fin n => i.val < w := by
    ext i
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    by_cases h : i.val < w
    · simp only [h, if_true, iff_true]
      intro hdvd
      have : (q : ℤ) ≤ 1 := Int.le_of_dvd one_pos hdvd
      have hq' : (2 : ℤ) ≤ (q : ℤ) := by exact_mod_cast hq
      linarith
    · simp only [h, if_false, iff_false, not_not]
      exact dvd_zero _
  rw [heq, card_lt_eq n w hw]

/-! ## Sufficiency -/

/-- **Sufficiency** (half of Proposition prop:spec-complete).
    For every q ≥ 2, every n, and every w ≤ n, there exists an integer
    lattice vector m ∈ ℤⁿ with |m|² ≤ n and parity weight w. -/
theorem completeness_sufficiency (q n w : ℕ) (hq : 2 ≤ q) (hw : w ≤ n) :
    ∃ m : Fin n → ℤ, normSq m ≤ (n : ℤ) ∧ parityWeight q m = w := by
  refine ⟨witnessVec n w, ?_, parityWeight_witness q n w hq hw⟩
  rw [normSq_witness n w hw]
  exact_mod_cast hw

/-! ## Tightness -/

/-- **Tightness** (half of Proposition prop:spec-complete).
    Any m ∈ ℤⁿ whose parity weight equals n has |m|² ≥ n.
    Equivalently: c_n(n-1) = 0 — no lattice vector of parity weight n
    can have squared norm strictly less than n. -/
theorem completeness_tightness (q n : ℕ) (_hq : 2 ≤ q)
    (m : Fin n → ℤ) (hpw : parityWeight q m = n) :
    (n : ℤ) ≤ normSq m := by
  -- Step 1: parity weight n means every coordinate is not divisible by q.
  have hfilt :
      (Finset.univ.filter fun i : Fin n => ¬ ((q : ℤ) ∣ m i))
        = (Finset.univ : Finset (Fin n)) := by
    apply Finset.eq_of_subset_of_card_le (Finset.filter_subset _ _) ?_
    rw [Finset.card_univ, Fintype.card_fin]
    have : (Finset.univ.filter fun i : Fin n => ¬ ((q : ℤ) ∣ m i)).card = n := hpw
    omega
  -- Filter = univ ⇒ each coordinate is nonzero (since q ∣ 0 is always true).
  have hall : ∀ i : Fin n, m i ≠ 0 := by
    intro i hi
    have hi_in : i ∈ Finset.univ.filter fun j : Fin n => ¬ ((q : ℤ) ∣ m j) := by
      rw [hfilt]; exact Finset.mem_univ i
    rw [Finset.mem_filter] at hi_in
    exact hi_in.2 (hi ▸ dvd_zero _)
  -- Step 2: each (m i)² ≥ 1 since m i is a nonzero integer.
  have hsq : ∀ i : Fin n, (1 : ℤ) ≤ (m i) ^ 2 := by
    intro i
    have hne := hall i
    rcases lt_trichotomy (m i) 0 with hlt | heq | hgt
    · nlinarith
    · exact (hne heq).elim
    · nlinarith
  -- Step 3: sum the bound.
  unfold normSq
  have hsum : ∑ _i : Fin n, (1 : ℤ) ≤ ∑ i, (m i) ^ 2 :=
    Finset.sum_le_sum (fun i _ => hsq i)
  have hconst : ∑ _i : Fin n, (1 : ℤ) = (n : ℤ) := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    simp
  linarith [hsum, hconst]

/-! ## Combined statement (Tier 1 main theorem) -/

/-- **Spectral completeness threshold (Proposition prop:spec-complete).**

    For every q ≥ 2 and every n ≥ 1, the threshold K = n is exactly
    where every parity weight class becomes nonempty:

    * Sufficiency: every w ∈ {0,…,n} is realised by some lattice
      vector with |m|² ≤ n.
    * Tightness:   no lattice vector of parity weight n has |m|² ≤ n - 1
      (equivalently, parity weight n forces |m|² ≥ n). -/
theorem spec_complete_threshold (q n : ℕ) (hq : 2 ≤ q) (_hn : 1 ≤ n) :
    (∀ w : ℕ, w ≤ n →
        ∃ m : Fin n → ℤ, normSq m ≤ (n : ℤ) ∧ parityWeight q m = w) ∧
    (∀ m : Fin n → ℤ, parityWeight q m = n → (n : ℤ) ≤ normSq m) :=
  ⟨fun w hw => completeness_sufficiency q n w hq hw,
   fun m hpw => completeness_tightness q n hq m hpw⟩

/-! ## Anchors at small (q, n) -/

/-- Anchor at q = 2, n = 4: the four canonical witnesses all have
    squared norm ≤ 4 and the expected parity weights. This pins the
    abstract sufficiency theorem to concrete numerical instances. -/
theorem spec_complete_n4_q2_anchor :
    normSq (witnessVec 4 0) = 0 ∧
    normSq (witnessVec 4 1) = 1 ∧
    normSq (witnessVec 4 2) = 2 ∧
    normSq (witnessVec 4 3) = 3 ∧
    normSq (witnessVec 4 4) = 4 ∧
    parityWeight 2 (witnessVec 4 0) = 0 ∧
    parityWeight 2 (witnessVec 4 1) = 1 ∧
    parityWeight 2 (witnessVec 4 2) = 2 ∧
    parityWeight 2 (witnessVec 4 3) = 3 ∧
    parityWeight 2 (witnessVec 4 4) = 4 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

end KstarFormal.SpecComplete
