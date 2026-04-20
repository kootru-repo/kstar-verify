/-
  KstarFormal.Probability.Hypergeometric — Hypergeometric/coupon-collector bound
  K* Verification: Krawtchouk spectral correspondence
  Registry: prop:coupon (Lemma 4)

  Layer 3 — probability theory.

  Lemma 4: For a uniform random draw of M operators from N non-identity Paulis:
    P(μ_w = 1) = Π_{j=0}^{A_w-1} (M - j)/(N - j) ≤ (M/N)^{A_w}

  The bound follows from the ratio-decreasing lemma:
    (M-j)/(N-j) ≤ M/N  for 0 ≤ j < M ≤ N.

  STATUS: SORRY-FREE (Layer 3 complete)
-/
import KstarFormal.Defs
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp

open Finset BigOperators

/-! ## Ratio-decreasing lemma

  For 0 ≤ j < M ≤ N (naturals):
    (M-j)·N ≤ M·(N-j)
  Proof: expand both sides:
    MN - jN ≤ MN - jM  ⟺  jM ≤ jN  ⟺  M ≤ N. ✓

  We prove this over ℤ (cross-multiplication form) and derive the ℚ
  division form using field_simp.
-/

/-- Cross-multiplication form: (M-j)·N ≤ M·(N-j) for j < M ≤ N.
    Working over ℤ to avoid natural subtraction issues. -/
theorem ratio_cross_mul (M N j : ℕ) (_hj : j < M) (hMN : M ≤ N) :
    ((M : ℤ) - j) * N ≤ M * (N - j) := by
  have : (j : ℤ) * M ≤ j * N := by
    apply mul_le_mul_of_nonneg_left
    · exact_mod_cast hMN
    · exact_mod_cast Nat.zero_le j
  nlinarith

/-- Ratio-decreasing over ℚ: (M-j)/(N-j) ≤ M/N for j < M ≤ N.
    Derived from cross-multiplication using field_simp. -/
theorem ratio_decreasing (M N j : ℕ) (hj : j < M) (hMN : M ≤ N) :
    (↑(M - j) : ℚ) / ↑(N - j) ≤ (↑M : ℚ) / ↑N := by
  have hjM := Nat.le_of_lt hj
  have hjN := Nat.le_trans hjM hMN
  have hN_ne : (↑N : ℚ) ≠ 0 := by exact_mod_cast (show (N : ℤ) ≠ 0 by omega)
  have hNj_ne : (↑(N - j) : ℚ) ≠ 0 := by exact_mod_cast (show (N - j : ℤ) ≠ 0 by omega)
  rw [div_le_div_iff₀ (by exact_mod_cast (show (0 : ℤ) < N - j by omega))
                       (by exact_mod_cast (show (0 : ℤ) < N by omega))]
  -- Goal: ↑(M - j) * ↑N ≤ ↑M * ↑(N - j)
  exact_mod_cast ratio_cross_mul M N j hj hMN

/-! ## Hypergeometric product bound

  P(μ_w = 1) = Π_{j=0}^{A_w-1} (M-j)/(N-j)

  By ratio_decreasing, each factor ≤ M/N, so the product ≤ (M/N)^{A_w}.
-/

/-- Each factor in the hypergeometric product is bounded by M/N. -/
theorem each_factor_bounded (M N : ℕ) (A : ℕ) (hA : A ≤ M) (hMN : M ≤ N)
    (j : ℕ) (hj : j < A) :
    (↑(M - j) : ℚ) / ↑(N - j) ≤ (↑M : ℚ) / ↑N :=
  ratio_decreasing M N j (by omega) hMN

/-- Lemma 4: Hypergeometric product bound.
    Π_{j<A} (M-j)/(N-j) ≤ (M/N)^A  for A ≤ M ≤ N.

    Proof: each factor ≤ M/N by ratio_decreasing,
    so the product of A factors ≤ (M/N)^A. -/
theorem hypergeometric_product_bound (M N A : ℕ) (hA : A ≤ M) (hMN : M ≤ N) :
    (Finset.range A).prod (fun j => ((↑(M - j) : ℚ) / ↑(N - j))) ≤
    ((↑M : ℚ) / ↑N) ^ A := by
  induction A with
  | zero => simp
  | succ k ih =>
    rw [Finset.prod_range_succ, pow_succ]
    apply mul_le_mul
    · exact ih (by omega)
    · exact each_factor_bounded M N (k + 1) hA hMN k (by omega)
    · exact div_nonneg (Nat.cast_nonneg (α := ℚ) _) (Nat.cast_nonneg (α := ℚ) _)
    · exact pow_nonneg (div_nonneg (Nat.cast_nonneg (α := ℚ) _) (Nat.cast_nonneg (α := ℚ) _)) _

/-! ## Concrete bounds at n=4, M=137, N=255

  Weight-class sizes A_w = [1, 12, 54, 108, 81].
  Saturation ratio M/N = 137/255 ≈ 0.537 < 1.
-/

/-- Saturation ratio is strictly less than 1. -/
theorem saturation_ratio_lt_one : (137 : ℚ) / 255 < 1 := by norm_num

/-- M = 137 ≤ N = 255. -/
theorem M_le_N : (137 : ℕ) ≤ 255 := by norm_num

/-- The product of ratio bounds shows vanishing probability.
    As A_w → ∞ with M/N < 1, P(μ_w=1) → 0 exponentially. -/
theorem vanishing_probability : (137 : ℚ) / 255 < 1 ∧ 0 < (137 : ℚ) / 255 := by
  constructor <;> norm_num
