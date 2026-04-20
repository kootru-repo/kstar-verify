/-
  KstarFormal.LinearAlgebra.Monotonicity — Eigenvalue monotonicity
  K* Verification: Krawtchouk spectral correspondence
  Registry: lem:monotone (Lemma 3)

  Lemma 3: λ_w(K) is non-decreasing in K for every weight class w.
  Concrete verification: compare eigenvalues at K=4 and K=5.

  Eigenvalues at K=4 (from c_w_K4 = [9, 8, 24, 32, 16]):
    λ = [144, 32, 64, 128, 256]
  Eigenvalues at K=5 (from c_w_K5 = [9, 56, 24, 32, 16]):
    λ = [144, 224, 64, 128, 256]

  Strict increase only at w=1 (shell k=5 contributes 48 vectors of weight 1).
-/
import KstarFormal.Combinatorics.GreedyRedist
import KstarFormal.LinearAlgebra.SpectralDecomp

/-! ## Eigenvalues at K=4 -/

def eigenvalues_K4 : List ℕ :=
  (List.range 5).map fun w =>
    gramEigenvalue_from_cw 4 (c_w_K4.getD w 0) w

theorem eigenvalues_K4_eq : eigenvalues_K4 = [144, 32, 64, 128, 256] := by
  native_decide

/-! ## Monotonicity: λ_w(K=4) ≤ λ_w(K=5) for all w -/

theorem eigenvalue_monotone_K4_K5 :
    ∀ w : Fin 5, eigenvalues_K4.getD w.val 0 ≤ eigenvalues_K5.getD w.val 0 := by
  intro w; fin_cases w <;> native_decide

/-! ## Strict increase at w=1

  The k=5 shell has 48 lattice vectors, all of parity weight 1.
  This is why c_1 jumps from 8 (at K=4) to 56 (at K=5):
    c_1(K=5) - c_1(K=4) = 56 - 8 = 48 = r_4(5) shell contribution to weight 1.
  Consequently, λ_1 jumps from 32 to 224.
-/

theorem c_w1_jump : c_w_K5.getD 1 0 - c_w_K4.getD 1 0 = 48 := by native_decide

theorem eigenvalue_w1_strict_increase :
    eigenvalues_K4.getD 1 0 < eigenvalues_K5.getD 1 0 := by native_decide

/-- No change at other weight classes (shell k=5 only contributes weight 1). -/
theorem eigenvalue_w0_unchanged :
    eigenvalues_K4.getD 0 0 = eigenvalues_K5.getD 0 0 := by native_decide
theorem eigenvalue_w2_unchanged :
    eigenvalues_K4.getD 2 0 = eigenvalues_K5.getD 2 0 := by native_decide
theorem eigenvalue_w3_unchanged :
    eigenvalues_K4.getD 3 0 = eigenvalues_K5.getD 3 0 := by native_decide
theorem eigenvalue_w4_unchanged :
    eigenvalues_K4.getD 4 0 = eigenvalues_K5.getD 4 0 := by native_decide

/-! ## Parity weight of shell k=5

  Shell k=5 contains 48 vectors. All have parity weight 1.
  Verification: c_w(K=5) - c_w(K=4) = [0, 48, 0, 0, 0].
-/

theorem shell5_parity_weight :
    (List.range 5).map (fun w => c_w_K5.getD w 0 - c_w_K4.getD w 0) = [0, 48, 0, 0, 0] := by
  native_decide
