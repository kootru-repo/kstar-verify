/-
  KstarFormal.Combinatorics.LatticeCountGeneric — Generic lattice-point counts
  K* Verification: Krawtchouk spectral correspondence

  Phase A (universal-n airtightness plan, 2026-04-17): generic parametric
  versions of
    * r_n(k)      — shell multiplicity, |{m ∈ ℤ^n : |m|² = k}|
    * N_n(K)      — ball count, |{m ∈ ℤ^n : |m|² ≤ K}|
    * c_w(K, n)   — parity-weight count in the ball
  for arbitrary dimension n, via structural recursion on n.

  The existing KstarFormal.Combinatorics.LatticeCount module defines these
  at n = 4 concrete with hardcoded tables. This module generalizes and
  proves n = 4 compatibility via native_decide.

  No new axioms; definitions are decidable by structural recursion.
-/
import KstarFormal.Combinatorics.LatticeCount
import Mathlib.Tactic.IntervalCases

/-! ## Generic shell count r_n(k)

  r_n(k) counts lattice vectors m ∈ ℤ^n with |m|² = k exactly.

  Recurrence (split on the last coordinate j):
    r_0(0)       = 1
    r_0(k+1)     = 0
    r_{n+1}(k)   = r_n(k) + 2 · Σ_{j ≥ 1, j² ≤ k} r_n(k - j²)

  The factor 2 captures the sign of j (±j); the j = 0 case is singular.
  Termination: structural recursion on the first (n) argument.
-/

/-- Generic shell count: |{m ∈ ℤ^n : |m|² = k}|. -/
def r_gen : ℕ → ℕ → ℕ
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k =>
      r_gen n k
        + 2 * ((Finset.range (Nat.sqrt k + 1)).filter
                 (fun j => 1 ≤ j ∧ j * j ≤ k)).sum
                 (fun j => r_gen n (k - j * j))

/-- Generic ball count: N_n(K) = Σ_{k = 0}^K r_n(k). -/
def ball_count_gen (n K : ℕ) : ℕ :=
  ((List.range (K + 1)).map (r_gen n)).sum

/-! ## n = 4 compatibility

  Every concrete lattice count used in Layer 1 is a specialization of the
  generic definition.  The specializations are decidable, so native_decide
  settles each within its finite range.
-/

/-- Generic shell count agrees with the concrete r4 on the range used by the paper. -/
theorem r_gen_eq_r4 : ∀ k, k ≤ 5 → r_gen 4 k = r4 k := by
  intro k hk
  interval_cases k <;> native_decide

/-- Generic ball count agrees with the concrete N4 on the range used by the paper. -/
theorem ball_count_gen_eq_N4 : ∀ K, K ≤ 5 → ball_count_gen 4 K = N4 K := by
  intro K hK
  interval_cases K <;> native_decide

/-- Concrete N4(5) = 137 via the generic route. -/
theorem ball_count_gen_4_5 : ball_count_gen 4 5 = 137 := by native_decide

/-! ## Generic parity-weight count c_w(K, n)

  c_w(K, n) counts lattice vectors m ∈ ℤ^n with
    |m|² ≤ K   and   (number of odd components of m) = w.

  Recurrence on n, splitting on the last coordinate j and tracking its
  parity:
    c_w(K, 0)       = [w = 0]
    c_w(K, n+1)     =   c_w(K, n)                               -- j = 0 (even)
                     + 2 · Σ_{j ≥ 1 even, j² ≤ K} c_w(K - j², n)
                     + 2 · Σ_{j ≥ 1 odd,  j² ≤ K} c_{w-1}(K - j², n)

  Termination: structural recursion on n (the first argument).
-/

/-- Generic parity-weight count. -/
def c_w_gen : ℕ → ℕ → ℕ → ℕ
  | 0, _K, 0 => 1
  | 0, _K, _ + 1 => 0
  | n + 1, K, w =>
      c_w_gen n K w
        + 2 * ((Finset.range (Nat.sqrt K + 1)).filter
                 (fun j => 1 ≤ j ∧ j * j ≤ K ∧ j % 2 = 0)).sum
                 (fun j => c_w_gen n (K - j * j) w)
        + 2 * (match w with
               | 0 => 0
               | wp + 1 =>
                   ((Finset.range (Nat.sqrt K + 1)).filter
                      (fun j => 1 ≤ j ∧ j * j ≤ K ∧ j % 2 = 1)).sum
                      (fun j => c_w_gen n (K - j * j) wp))

/-- At n = 4, K = 5, the generic c_w_gen reproduces the certified list. -/
theorem c_w_gen_4_5_eq_c_w_K5 :
    [c_w_gen 4 5 0, c_w_gen 4 5 1, c_w_gen 4 5 2, c_w_gen 4 5 3, c_w_gen 4 5 4]
      = c_w_K5 := by
  native_decide

/-- Every concrete parity-weight count at n = 4, K = 5 matches the certified list. -/
theorem c_w_gen_4_5_matches :
    c_w_gen 4 5 0 = 9 ∧ c_w_gen 4 5 1 = 56 ∧ c_w_gen 4 5 2 = 24
      ∧ c_w_gen 4 5 3 = 32 ∧ c_w_gen 4 5 4 = 16 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> native_decide

/-- Generic ball count decomposes into parity-weight-w sums (for any concrete n, K
    within the paper's regime; proved at n = 4, K = 5 via native_decide). -/
theorem c_w_gen_sum_eq_ball_count_n4_K5 :
    (c_w_gen 4 5 0) + (c_w_gen 4 5 1) + (c_w_gen 4 5 2)
      + (c_w_gen 4 5 3) + (c_w_gen 4 5 4) = ball_count_gen 4 5 := by
  native_decide

/-! ## Generic Krawtchouk eigenvalue

  The existing `gramEigenvalue_from_cw` in SpectralDecomp.lean is already
  generic in n.  Nothing new is needed here; we re-expose the generic
  eigenvalue computation with the generic c_w_gen input for clarity.
-/

/-- Generic Krawtchouk eigenvalue via the generic parity-weight count.
    λ_w(K; n) = 2^n · c_w(K, n) / C(n, w). -/
def krawtchouk_eigenvalue_gen (n K w : ℕ) : ℕ :=
  2 ^ n * c_w_gen n K w / Nat.choose n w

/-- At n = 4, K = 5, the generic eigenvalue formula reproduces the certified list. -/
theorem krawtchouk_eigenvalue_gen_4_5 :
    [krawtchouk_eigenvalue_gen 4 5 0, krawtchouk_eigenvalue_gen 4 5 1,
     krawtchouk_eigenvalue_gen 4 5 2, krawtchouk_eigenvalue_gen 4 5 3,
     krawtchouk_eigenvalue_gen 4 5 4] = [144, 224, 64, 128, 256] := by
  native_decide

/-- All Krawtchouk eigenvalues at n = 4, K = 5 are strictly positive. -/
theorem krawtchouk_eigenvalue_gen_4_5_pos :
    ∀ w : Fin 5, 0 < krawtchouk_eigenvalue_gen 4 5 w.val := by
  intro w; fin_cases w <;> native_decide
