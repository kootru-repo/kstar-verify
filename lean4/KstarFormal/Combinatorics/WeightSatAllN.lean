/-
  KstarFormal.Combinatorics.WeightSatAllN — Weight-1,2 saturation
  at K = 5 for all n ≥ 4 (Tier 1B).
  K* Verification: Krawtchouk spectral correspondence

  Formalizes the manuscript claim that K* = 5 saturates Pauli weight
  classes 1 and 2 for every n ≥ 4 after greedy redistribution of
  the parity-weight counts. The closed-form lattice counts at K = 5 are:

      c_0(5, n) = 1 + 2n
      c_1(5, n) = 4n² - 2n           (for n ≥ 1)
      c_2(5, n) = 2n (n - 1)

  These are derived by direct enumeration of m ∈ ℤⁿ with |m|² ≤ 5
  (each |m_i| ∈ {0, 1, 2}, with parity-weight w meaning exactly w
  coordinates are odd, i.e., ±1).

  The Pauli totals at q = 2 are A_w(n) = 3^w · C(n, w):
      A_0(n) = 1, A_1(n) = 3n, A_2(n) = 9·n·(n-1)/2.

  The load-bearing inequality, after the greedy redistribution that
  fills A_0 and A_1 from the surplus of c_0 and c_1, is that the
  remaining lattice mass available for w = 2 is at least A_2:

      c_0 + c_1 + c_2 - A_0 - A_1  ≥  A_2.

  We prove this with `nlinarith` for all n ≥ 1, which a fortiori
  gives the manuscript claim for n ≥ 4.

  Status: Tier 1B formalization, sorry-free.
-/
import KstarFormal.Combinatorics.GreedyRedist
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Zify

namespace KstarFormal.WeightSatAllN

/-! ## Closed-form lattice counts at K = 5 -/

/-- Number of m ∈ ℤⁿ with |m|² ≤ 5 and parity weight 0
    (all coordinates even, i.e., in {0, ±2}, sum-of-squares ≤ 5).
    Equals 1 + 2n: the all-zero vector plus, for each of n positions,
    a single ±2 in that position. -/
def c0_K5 (n : ℕ) : ℕ := 1 + 2 * n

/-- Number of m ∈ ℤⁿ with |m|² ≤ 5 and parity weight 1.
    Pick the odd-coord position (n choices, sign ±1: 2 choices),
    then n-1 even-coord choices with sum-of-squares ≤ 4
    (1 all-zero + 2(n-1) configurations with one ±2).
    Total: 2n · (1 + 2(n-1)) = 2n · (2n-1) = 4n² - 2n. -/
def c1_K5 (n : ℕ) : ℕ := 4 * n * n - 2 * n

/-- Number of m ∈ ℤⁿ with |m|² ≤ 5 and parity weight 2.
    Pick 2 odd-coord positions (C(n,2) choices, signs ±1: 4 choices),
    then n-2 even-coord choices with sum-of-squares ≤ 3
    (only the all-zero configuration: ±2 contributes 4 > 3).
    Total: 4 · C(n,2) = 2n(n-1). -/
def c2_K5 (n : ℕ) : ℕ := 2 * n * (n - 1)

/-! ## Pauli weight-class totals A_w(n) at q = 2 -/

/-- A_0(n) = 1: the identity operator. -/
def A0_n : ℕ := 1

/-- A_1(n) = 3n: three Pauli types per single qubit. -/
def A1_n (n : ℕ) : ℕ := 3 * n

/-- 2·A_2(n) = 9·n·(n-1).  We work with twice A_2 to avoid the
    natural-number division in 9·C(n,2) = 9n(n-1)/2. -/
def twoA2_n (n : ℕ) : ℕ := 9 * n * (n - 1)

/-! ## Anchor at n = 4

  The closed forms must reproduce the certified `c_w_K5` list at n = 4
  (verified elsewhere via `native_decide` against direct lattice enumeration).
-/

/-- At n = 4, the closed-form counts agree with the certified list. -/
theorem c_K5_n4_anchor :
    c0_K5 4 = 9 ∧ c1_K5 4 = 56 ∧ c2_K5 4 = 24 := by
  refine ⟨?_, ?_, ?_⟩
  · rfl
  · rfl
  · rfl

/-- At n = 4, the Pauli totals match the existing `A_w_n4` list. -/
theorem A_n4_anchor :
    A0_n = 1 ∧ A1_n 4 = 12 ∧ twoA2_n 4 = 108 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-! ## The load-bearing inequality

  The greedy redistribution starting from `(c_0, c_1, c_2)` and
  filling `(A_0, A_1, A_2)` in order succeeds iff the cumulative
  lattice mass at each prefix is at least the cumulative budget at
  that prefix:

      c_0           ≥ A_0
      c_0 + c_1     ≥ A_0 + A_1
      c_0 + c_1 + c_2 ≥ A_0 + A_1 + A_2.

  At K = 5 the first two are immediate (c_0 = 1 + 2n ≥ 1 = A_0;
  c_0 + c_1 = 1 + 2n + 4n² - 2n = 1 + 4n² ≥ 1 + 3n = A_0 + A_1
  for n ≥ 1).  The third is the load-bearing claim, equivalent
  after multiplying by 2 (to clear the half in A_2) to

      2(1 + 2n + 4n² - 2n + 2n² - 2n) ≥ 2 + 6n + 9n(n-1)
      ⟺ 12n² - 4n + 2                 ≥ 9n² - 3n + 2
      ⟺ 3n² ≥ n,

  which holds for every n ≥ 1.
-/

/-- Step 1: c_0 ≥ A_0 (trivially, 1 + 2n ≥ 1). -/
theorem step0_K5 (n : ℕ) : c0_K5 n ≥ A0_n := by
  unfold c0_K5 A0_n
  omega

/-- Step 2: c_0 + c_1 ≥ A_0 + A_1, i.e. 1 + 4n² ≥ 1 + 3n,
    valid for all n (and a fortiori for n ≥ 1). -/
theorem step1_K5 (n : ℕ) (hn : 1 ≤ n) : c0_K5 n + c1_K5 n ≥ A0_n + A1_n n := by
  unfold c0_K5 c1_K5 A0_n A1_n
  have h_sub : 2 * n ≤ 4 * n * n := by nlinarith
  zify [h_sub]
  nlinarith [sq_nonneg ((n : ℤ) - 1), hn]

/-- Step 3 (the load-bearing inequality): the doubled cumulative mass
    of c_0, c_1, c_2 dominates the doubled cumulative budget A_0, A_1, A_2.

    After multiplying through to avoid the division in A_2 = 9n(n-1)/2,
    the inequality reduces to 3n² ≥ n. -/
theorem step2_K5 (n : ℕ) (hn : 1 ≤ n) :
    2 * (c0_K5 n + c1_K5 n + c2_K5 n) ≥ 2 * A0_n + 2 * A1_n n + twoA2_n n := by
  unfold c0_K5 c1_K5 c2_K5 A0_n A1_n twoA2_n
  -- Reduce nat subtraction by exposing the underlying inequalities.
  -- 4n² ≥ 2n (from n ≥ 1, used inside c1_K5) and (n-1) handled by omega.
  -- Eliminate nat subtraction in c1, then cast to ℤ to handle (n-1).
  have h1 : 2 * n ≤ 4 * n * n := by nlinarith
  zify [h1, hn]
  nlinarith [sq_nonneg ((n : ℤ) - 1), hn, sq_nonneg (n : ℤ)]

/-! ## Tier 1B main theorem -/

/-- **Tier 1B main theorem.** For all n ≥ 1 (a fortiori for n ≥ 4),
    the K = 5 lattice provides enough cumulative mass through parity
    weights 0, 1, 2 to saturate the Pauli weight totals A_0, A_1, A_2
    via greedy redistribution. Equivalently: the doubled cumulative
    surplus is non-negative at every step.

    This is the load-bearing fact behind the manuscript claim
    "K* = 5 saturates weight classes 1 and 2 for every n ≥ 4". -/
theorem weight12_saturation_K5_all_n (n : ℕ) (hn : 1 ≤ n) :
    c0_K5 n ≥ A0_n ∧
    c0_K5 n + c1_K5 n ≥ A0_n + A1_n n ∧
    2 * (c0_K5 n + c1_K5 n + c2_K5 n) ≥ 2 * A0_n + 2 * A1_n n + twoA2_n n :=
  ⟨step0_K5 n, step1_K5 n hn, step2_K5 n hn⟩

/-- Specialization for the manuscript-stated regime n ≥ 4. -/
theorem weight12_saturation_K5_nge4 (n : ℕ) (hn : 4 ≤ n) :
    c0_K5 n ≥ A0_n ∧
    c0_K5 n + c1_K5 n ≥ A0_n + A1_n n ∧
    2 * (c0_K5 n + c1_K5 n + c2_K5 n) ≥ 2 * A0_n + 2 * A1_n n + twoA2_n n :=
  weight12_saturation_K5_all_n n (by omega)

/-- Cross-check at the n = 4 anchor: the existing certified
    `M_w_K5 = [1, 12, 54, 54, 16]` saturates A_1 and A_2. This
    pins the abstract Tier 1B claim to the concrete Lean infrastructure. -/
theorem n4_anchor_consistency :
    M_w_K5.getD 1 0 = A_w_n4.getD 1 0 ∧
    M_w_K5.getD 2 0 = A_w_n4.getD 2 0 := by
  refine ⟨?_, ?_⟩
  · native_decide
  · native_decide

/-! ## Multi-n specialization sanity checks (technique #4)

  Concrete numerical checks of the closed forms `c0_K5`, `c1_K5`, `c2_K5`,
  `A1_n`, `twoA2_n` at n ∈ {4, 5, 6, 7, 10}. These pin the abstract
  symbolic identities to literal values that can be cross-validated against
  the Python harness (`cross_validate_lean.py`, Section 1) and against any
  textbook lattice-point computation. If a refactor silently changes the
  algebra, these `decide`-based anchors fail loudly.
-/

-- n = 5
theorem cK5_n5  : c0_K5 5  =  11 ∧ c1_K5 5  =  90 ∧ c2_K5 5  =  40 := by decide
theorem AK5_n5  : A1_n 5   =  15 ∧ twoA2_n 5 = 180 := by decide
-- n = 6
theorem cK5_n6  : c0_K5 6  =  13 ∧ c1_K5 6  = 132 ∧ c2_K5 6  =  60 := by decide
theorem AK5_n6  : A1_n 6   =  18 ∧ twoA2_n 6 = 270 := by decide
-- n = 7
theorem cK5_n7  : c0_K5 7  =  15 ∧ c1_K5 7  = 182 ∧ c2_K5 7  =  84 := by decide
theorem AK5_n7  : A1_n 7   =  21 ∧ twoA2_n 7 = 378 := by decide
-- n = 10
theorem cK5_n10 : c0_K5 10 =  21 ∧ c1_K5 10 = 380 ∧ c2_K5 10 = 180 := by decide
theorem AK5_n10 : A1_n 10  =  30 ∧ twoA2_n 10 = 810 := by decide

-- Note: a former `step2_K5_concrete` inequality-form anchor was removed
-- after constant-fuzzing showed its slack let nearby literal perturbations
-- still close via `decide`. The strict equality form `step2_K5_surplus`
-- below pins the surplus exactly to `3n^2 - n` and is the authoritative
-- numerical anchor.

/-- Negative anchor: at the *boundary* the inequality must be tight or
    strict by a positive surplus, never reversed. We pin the surplus
    `2(c_0+c_1+c_2) - (2A_0+2A_1+2A_2)` at the spot-checked n's. The
    closed-form value is `3n² - n`, which we verify numerically. Any
    sign flip in the algebra would fail one of these. -/
theorem step2_K5_surplus :
    2 * (c0_K5 5  + c1_K5 5  + c2_K5 5)  - (2 * A0_n + 2 * A1_n 5  + twoA2_n 5)  =  70  ∧
    2 * (c0_K5 6  + c1_K5 6  + c2_K5 6)  - (2 * A0_n + 2 * A1_n 6  + twoA2_n 6)  = 102  ∧
    2 * (c0_K5 7  + c1_K5 7  + c2_K5 7)  - (2 * A0_n + 2 * A1_n 7  + twoA2_n 7)  = 140  ∧
    2 * (c0_K5 10 + c1_K5 10 + c2_K5 10) - (2 * A0_n + 2 * A1_n 10 + twoA2_n 10) = 290 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-! ## Weight-3 closed form and strict undersaturation

  Extending the analysis to parity weight 3, the closed-form lattice
  count at K = 5 is

      c_3(5, n) = 8 · C(n, 3)

  derived by picking 3 odd-coord positions (`C(n,3)` choices), assigning
  ±1 to each (2³ = 8 sign choices), and observing that the remaining
  n - 3 even coordinates must all be zero because any nonzero ±2
  contributes 4 > 5 - 3 to the squared norm.

  The Pauli weight-3 budget is `A_3(n) = 27 · C(n, 3)`. The closed-form
  identity `c_3(5, n) + 19 · C(n, 3) = A_3(n)` exhibits the deficit
  exactly: 19 · C(n, 3) operators per n. For every n ≥ 3 the deficit is
  strictly positive, so **K = 5 never saturates weight 3 for n ≥ 3**.
  This is the load-bearing fact behind the manuscript-stated scope of
  the K* = 5 saturation result (weights 1 and 2 only). -/

/-- c_3(5, n) = 8 · C(n, 3): pick 3 odd positions, ±1 each, even
    coordinates forced to zero (any ±2 contributes 4 > 2 = 5 - 3). -/
def c3_K5 (n : ℕ) : ℕ := 8 * Nat.choose n 3

/-- A_3(n) = 27 · C(n, 3): three Pauli types (X, Y, Z) per qubit on a
    chosen 3-element subset. -/
def A3_n (n : ℕ) : ℕ := 27 * Nat.choose n 3

/-- Anchor at n = 4: c_3(5, 4) = 32 (matches certified `c_w_K5[3] = 32`). -/
theorem c3_K5_n4 : c3_K5 4 = 32 := by decide

/-- Anchor at n = 4: A_3(4) = 108 (= 27 · C(4,3) = 27 · 4). -/
theorem A3_n4 : A3_n 4 = 108 := by decide

/-- Multi-n closed-form spot checks (techniques #4 and constant fuzzing). -/
theorem c3_K5_spot :
    c3_K5 3 = 8   ∧ c3_K5 4 = 32  ∧ c3_K5 5 = 80  ∧
    c3_K5 6 = 160 ∧ c3_K5 7 = 280 ∧ c3_K5 10 = 960 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- Multi-n A_3 spot checks. -/
theorem A3_n_spot :
    A3_n 3 = 27   ∧ A3_n 4 = 108  ∧ A3_n 5 = 270  ∧
    A3_n 6 = 540  ∧ A3_n 7 = 945  ∧ A3_n 10 = 3240 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- **Closed-form deficit identity.** The K = 5 weight-3 lattice count
    differs from the Pauli weight-3 budget by exactly `19 · C(n, 3)` for
    every n. -/
theorem weight3_deficit_identity (n : ℕ) :
    c3_K5 n + 19 * Nat.choose n 3 = A3_n n := by
  unfold c3_K5 A3_n
  omega

/-- **K = 5 strictly under-saturates weight 3 for all n ≥ 3.**
    This is the closed-form extension of the saturation result: where
    weights 1 and 2 saturate for every n ≥ 1, weight 3 fails by
    `19 · C(n, 3) > 0` operators for every n ≥ 3. -/
theorem weight3_strict_undersaturation_K5 (n : ℕ) (hn : 3 ≤ n) :
    c3_K5 n < A3_n n := by
  unfold c3_K5 A3_n
  have hpos : 0 < Nat.choose n 3 := Nat.choose_pos hn
  omega

end KstarFormal.WeightSatAllN