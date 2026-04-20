/-
  KstarFormal.Combinatorics.DelsarteDual — Delsarte LP dual feasibility certificate at K*=5
  K* Verification: Krawtchouk spectral correspondence

  In the Delsarte LP framework for association schemes, the Gram matrix
  G(K) lies in the Bose-Mesner algebra of H(n,q). Its spectral decomposition
  G(K) = Σ_w λ_w E_w has eigenvalues λ_w = d · f_w where f_w = c_w / C(n,w).

  The LP dual certificate for G(K) ≥ 0 (PSD) is: f_w ≥ 0 for all w.
  The LP dual certificate for full rank is: f_w > 0 for all w.

  This file computes the normalized dual variables f_w over Q and verifies:
    1. All f_w > 0 (dual strict feasibility, i.e., full rank)
    2. The dual objective Σ f_w C(n,w) = N_4(5) = 137
    3. Eigenvalue recovery λ_w = d · f_w
    4. The Krawtchouk spectral coefficients partition correctly

  All proofs by native_decide over Q.

  Status: Tier 1 formalization, sorry-free.
-/
import KstarFormal.Defs
import KstarFormal.Combinatorics.LatticeCount
import KstarFormal.LinearAlgebra.SpectralDecomp

open Finset BigOperators

/-! ## Delsarte dual variables (normalized spectral coefficients)

  f_w = c_w(K) / C(n, w) is the "per-weight-class" lattice density.
  For n=4, K=5: f_w = [9/1, 56/4, 24/6, 32/4, 16/1] = [9, 14, 4, 8, 16].
  All integer-valued (a consequence of the divisibility structure of Z^n lattices).
-/

/-- Delsarte dual variables at K=5: f_w = c_w / C(4,w).
    These are the normalized spectral coefficients of the Gram matrix. -/
def delsarteDual_K5 : List ℕ :=
  (List.range 5).map fun w =>
    c_w_K5.getD w 0 / Nat.choose 4 w

theorem delsarteDual_K5_eq : delsarteDual_K5 = [9, 14, 4, 8, 16] := by native_decide

/-! ## Dual strict feasibility: all f_w > 0 -/

theorem delsarteDual_K5_all_pos :
    ∀ w : Fin 5, delsarteDual_K5.getD w.val 0 > 0 := by
  intro w; fin_cases w <;> native_decide

/-- No zero dual variables (contrasts with K=3 where f_4 = 0). -/
theorem delsarteDual_K5_no_zero :
    (List.range 5).Forall fun w => delsarteDual_K5.getD w 0 ≠ 0 := by native_decide

/-! ## Delsarte dual variables at K=3 (rank-deficient case) -/

/-- Parity weight counts at K=3 (local copy; canonical in FidelityDichotomy.lean). -/
private def c_w_K3_local : List ℕ := [1, 8, 24, 32, 0]

/-- At K=3, f_4 = c_4/C(4,4) = 0/1 = 0, so dual is not strictly feasible. -/
def delsarteDual_K3 : List ℕ :=
  (List.range 5).map fun w =>
    c_w_K3_local.getD w 0 / Nat.choose 4 w

theorem delsarteDual_K3_eq : delsarteDual_K3 = [1, 2, 4, 8, 0] := by native_decide

/-- K=3 fails strict feasibility: f_4 = 0 (the Delsarte obstruction). -/
theorem delsarteDual_K3_has_zero :
    delsarteDual_K3.getD 4 0 = 0 := by native_decide

/-- The Delsarte obstruction at K=3 is UNIQUE to weight 4. -/
theorem delsarteDual_K3_unique_zero :
    (∀ w : Fin 4, delsarteDual_K3.getD w.val 0 > 0) ∧
    delsarteDual_K3.getD 4 0 = 0 := by
  constructor
  · intro w; fin_cases w <;> native_decide
  · native_decide

/-! ## Dual objective: Σ f_w · C(n,w) = N_4(K)

  The LP dual objective recovers the total measurement count.
  This is the Delsarte LP bound: |S| = Σ_w f_w · C(n,w).
-/

theorem delsarteDual_objective_K5 :
    ((List.range 5).map fun w =>
      delsarteDual_K5.getD w 0 * Nat.choose 4 w).sum = 137 := by native_decide

theorem delsarteDual_objective_K3 :
    ((List.range 5).map fun w =>
      delsarteDual_K3.getD w 0 * Nat.choose 4 w).sum = 65 := by native_decide

/-! ## Eigenvalue recovery: λ_w = d · f_w

  The Gram eigenvalues are d times the dual variables.
  This connects the LP dual certificate to the spectral decomposition.
-/

theorem delsarte_eigenvalue_recovery_K5 :
    (List.range 5).map (fun w => 16 * delsarteDual_K5.getD w 0)
    = eigenvalues_K5 := by native_decide

/-- Alternative verification: eigenvalue formula is consistent.
    λ_w = 2^n · c_w / C(n,w) = d · f_w. -/
theorem delsarte_eigenvalue_consistent :
    ∀ w : Fin 5,
      gramEigenvalue_from_cw 4 (c_w_K5.getD w.val 0) w.val
      = 16 * delsarteDual_K5.getD w.val 0 := by
  intro w; fin_cases w <;> native_decide

/-! ## Dual monotonicity: K=3 → K=5

  The Delsarte dual variables grow monotonically with K.
  f_w(K=3) ≤ f_w(K=5) for all w.
-/

theorem delsarteDual_monotone_K3_K5 :
    ∀ w : Fin 5, delsarteDual_K3.getD w.val 0 ≤ delsarteDual_K5.getD w.val 0 := by
  intro w; fin_cases w <;> native_decide

/-! ## Delsarte LP optimality gap

  The difference between K=4 (first full rank) and K*=5 (dynamic threshold)
  in terms of dual variables quantifies the "margin" K* provides.
-/

/-- Parity weight counts at K=4 (local copy; canonical in GreedyRedist.lean). -/
private def c_w_K4_local : List ℕ := [9, 8, 24, 32, 16]

/-- Delsarte dual at K=4. -/
def delsarteDual_K4 : List ℕ :=
  (List.range 5).map fun w =>
    c_w_K4_local.getD w 0 / Nat.choose 4 w

theorem delsarteDual_K4_eq : delsarteDual_K4 = [9, 2, 4, 8, 16] := by native_decide

/-- K=4 also has strictly positive duals (full rank starts at K=4). -/
theorem delsarteDual_K4_all_pos :
    ∀ w : Fin 5, delsarteDual_K4.getD w.val 0 > 0 := by
  intro w; fin_cases w <;> native_decide

/-- The K=4 → K=5 improvement is concentrated in weight 1.
    f_1(K=4)=2, f_1(K=5)=14: 7x increase from the 48 shell-5 vectors. -/
theorem delsarte_weight1_jump :
    delsarteDual_K5.getD 1 0 - delsarteDual_K4.getD 1 0 = 12 ∧
    delsarteDual_K5.getD 0 0 = delsarteDual_K4.getD 0 0 ∧
    delsarteDual_K5.getD 2 0 = delsarteDual_K4.getD 2 0 ∧
    delsarteDual_K5.getD 3 0 = delsarteDual_K4.getD 3 0 ∧
    delsarteDual_K5.getD 4 0 = delsarteDual_K4.getD 4 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> native_decide
