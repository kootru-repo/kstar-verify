/-
  KstarFormal.Combinatorics.KstarLabels — Explicit 137-label enumeration
  K* Verification: Krawtchouk spectral correspondence

  Deterministic extraction of all 137 lattice points m ∈ Z⁴ with |m|² ≤ 5.
  Verifies:
    1. Exactly 137 points (matching N_4(5))
    2. Parity-weight distribution c_w = [9, 56, 24, 32, 16]
    3. Each shell count matches r_4(k)

  All proofs by native_decide over the finite enumeration.

  Status: Tier 1D formalization, sorry-free.
-/
import Mathlib.Data.Int.Interval
import KstarFormal.Combinatorics.LatticeCount

open Finset

/-! ## Lattice point enumeration in Z⁴ -/

/-- Component range: integers from -2 to 2 (sufficient since m_i² ≤ 5 implies |m_i| ≤ 2). -/
private abbrev R : Finset Int := Finset.Icc (-2) 2

/-- Squared norm of a 4-tuple. -/
private def sqNorm (m : ((Int × Int) × Int) × Int) : Int :=
  let ⟨⟨⟨a, b⟩, c⟩, d⟩ := m
  a * a + b * b + c * c + d * d

/-- Parity weight: number of odd components. -/
def parityWeight4 (m : ((Int × Int) × Int) × Int) : Nat :=
  let ⟨⟨⟨a, b⟩, c⟩, d⟩ := m
  (if a % 2 ≠ 0 then 1 else 0) + (if b % 2 ≠ 0 then 1 else 0) +
  (if c % 2 ≠ 0 then 1 else 0) + (if d % 2 ≠ 0 then 1 else 0)

/-- All lattice points m ∈ Z⁴ with |m|² ≤ 5 (the K*=5 measurement labels). -/
def kstarLabels_n4 : Finset (((Int × Int) × Int) × Int) :=
  (((R ×ˢ R) ×ˢ R) ×ˢ R).filter fun m => decide (sqNorm m ≤ 5)

/-! ## Cardinality: exactly 137 labels -/

theorem kstarLabels_card : kstarLabels_n4.card = 137 := by native_decide

/-! ## Parity-weight distribution matches c_w(K=5) = [9, 56, 24, 32, 16] -/

theorem kstarLabels_pw0 :
    (kstarLabels_n4.filter fun m => parityWeight4 m = 0).card = 9 := by native_decide

theorem kstarLabels_pw1 :
    (kstarLabels_n4.filter fun m => parityWeight4 m = 1).card = 56 := by native_decide

theorem kstarLabels_pw2 :
    (kstarLabels_n4.filter fun m => parityWeight4 m = 2).card = 24 := by native_decide

theorem kstarLabels_pw3 :
    (kstarLabels_n4.filter fun m => parityWeight4 m = 3).card = 32 := by native_decide

theorem kstarLabels_pw4 :
    (kstarLabels_n4.filter fun m => parityWeight4 m = 4).card = 16 := by native_decide

/-- Full distribution matches the certified c_w_K5 vector. -/
theorem kstarLabels_cw_match :
    [(kstarLabels_n4.filter fun m => parityWeight4 m = 0).card,
     (kstarLabels_n4.filter fun m => parityWeight4 m = 1).card,
     (kstarLabels_n4.filter fun m => parityWeight4 m = 2).card,
     (kstarLabels_n4.filter fun m => parityWeight4 m = 3).card,
     (kstarLabels_n4.filter fun m => parityWeight4 m = 4).card] = c_w_K5 := by native_decide

/-! ## Shell decomposition: each shell has the correct count -/

/-- Points on shell k (|m|² = k). -/
def shellPoints (k : Int) : Finset (((Int × Int) × Int) × Int) :=
  kstarLabels_n4.filter fun m => sqNorm m = k

theorem shell0_card : (shellPoints 0).card = 1 := by native_decide
theorem shell1_card : (shellPoints 1).card = 8 := by native_decide
theorem shell2_card : (shellPoints 2).card = 24 := by native_decide
theorem shell3_card : (shellPoints 3).card = 32 := by native_decide
theorem shell4_card : (shellPoints 4).card = 24 := by native_decide
theorem shell5_card : (shellPoints 5).card = 48 := by native_decide

/-- Shell counts match r_4(k). -/
theorem shells_match_r4 :
    [(shellPoints 0).card, (shellPoints 1).card, (shellPoints 2).card,
     (shellPoints 3).card, (shellPoints 4).card, (shellPoints 5).card]
    = [r4 0, r4 1, r4 2, r4 3, r4 4, r4 5] := by native_decide
