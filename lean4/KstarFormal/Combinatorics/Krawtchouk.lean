/-
  KstarFormal.Combinatorics.Krawtchouk — Krawtchouk polynomial definition and properties
  K* Verification: Krawtchouk spectral correspondence
  Registry: prop:spectral_q_main (Lemma 5), foundation for all spectral results
-/
import Mathlib.Data.Int.Basic
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.FinCases


open Finset BigOperators

/-! ## Definition -/

/-- Binary Krawtchouk polynomial K_w(h; n) evaluated at concrete values.
    Uses an explicit recursive definition for decidability. -/
def krawtchouk_term (n w h j : ℕ) : ℤ :=
  (-1 : ℤ) ^ j * (Nat.choose h j : ℤ) * (Nat.choose (n - h) (w - j) : ℤ)

/-- Binary Krawtchouk polynomial K_w(h; n) for small w (unrolled). -/
def krawtchouk (n w h : ℕ) : ℤ :=
  (List.range (w + 1)).map (krawtchouk_term n w h) |>.sum

/-! ## Concrete values for n=4 -/

/-- Krawtchouk matrix entry K_w(h; 4) as a function on Fin 5 × Fin 5. -/
def krawtchouk_n4 (w h : Fin 5) : ℤ := krawtchouk 4 w.val h.val

-- Verify all 25 entries of the 5×5 Krawtchouk matrix via native_decide.
-- Row 0: [1, 1, 1, 1, 1]
-- Row 1: [4, 2, 0, -2, -4]
-- Row 2: [6, 0, -2, 0, 6]
-- Row 3: [4, -2, 0, 2, -4]
-- Row 4: [1, -1, 1, -1, 1]

theorem kraw_0_0 : krawtchouk 4 0 0 = 1 := by native_decide
theorem kraw_0_1 : krawtchouk 4 0 1 = 1 := by native_decide
theorem kraw_0_2 : krawtchouk 4 0 2 = 1 := by native_decide
theorem kraw_0_3 : krawtchouk 4 0 3 = 1 := by native_decide
theorem kraw_0_4 : krawtchouk 4 0 4 = 1 := by native_decide

theorem kraw_1_0 : krawtchouk 4 1 0 = 4 := by native_decide
theorem kraw_1_1 : krawtchouk 4 1 1 = 2 := by native_decide
theorem kraw_1_2 : krawtchouk 4 1 2 = 0 := by native_decide
theorem kraw_1_3 : krawtchouk 4 1 3 = -2 := by native_decide
theorem kraw_1_4 : krawtchouk 4 1 4 = -4 := by native_decide

theorem kraw_2_0 : krawtchouk 4 2 0 = 6 := by native_decide
theorem kraw_2_1 : krawtchouk 4 2 1 = 0 := by native_decide
theorem kraw_2_2 : krawtchouk 4 2 2 = -2 := by native_decide
theorem kraw_2_3 : krawtchouk 4 2 3 = 0 := by native_decide
theorem kraw_2_4 : krawtchouk 4 2 4 = 6 := by native_decide

theorem kraw_3_0 : krawtchouk 4 3 0 = 4 := by native_decide
theorem kraw_3_1 : krawtchouk 4 3 1 = -2 := by native_decide
theorem kraw_3_2 : krawtchouk 4 3 2 = 0 := by native_decide
theorem kraw_3_3 : krawtchouk 4 3 3 = 2 := by native_decide
theorem kraw_3_4 : krawtchouk 4 3 4 = -4 := by native_decide

theorem kraw_4_0 : krawtchouk 4 4 0 = 1 := by native_decide
theorem kraw_4_1 : krawtchouk 4 4 1 = -1 := by native_decide
theorem kraw_4_2 : krawtchouk 4 4 2 = 1 := by native_decide
theorem kraw_4_3 : krawtchouk 4 4 3 = -1 := by native_decide
theorem kraw_4_4 : krawtchouk 4 4 4 = 1 := by native_decide

/-! ## Orthogonality verification

  Krawtchouk orthogonality (binary case):
    Σ_{h=0}^{4} C(4,h) · K_w(h;4) · K_{w'}(h;4) = 16 · C(4,w) · δ_{w,w'}

  We verify all 25 (w,w') pairs. Diagonal entries should be 16·C(4,w).
  Off-diagonal entries should be 0.
-/

/-- Weighted inner product of Krawtchouk rows w and w'. -/
def krawtchouk_inner (w w' : ℕ) : ℤ :=
  (List.range 5).map (fun h =>
    (Nat.choose 4 h : ℤ) * krawtchouk 4 w h * krawtchouk 4 w' h) |>.sum

-- Diagonal: ⟨K_w, K_w⟩ = 16 · C(4,w)
theorem ortho_diag_0 : krawtchouk_inner 0 0 = 16 * 1 := by native_decide
theorem ortho_diag_1 : krawtchouk_inner 1 1 = 16 * 4 := by native_decide
theorem ortho_diag_2 : krawtchouk_inner 2 2 = 16 * 6 := by native_decide
theorem ortho_diag_3 : krawtchouk_inner 3 3 = 16 * 4 := by native_decide
theorem ortho_diag_4 : krawtchouk_inner 4 4 = 16 * 1 := by native_decide

-- Off-diagonal: ⟨K_w, K_{w'}⟩ = 0
theorem ortho_01 : krawtchouk_inner 0 1 = 0 := by native_decide
theorem ortho_02 : krawtchouk_inner 0 2 = 0 := by native_decide
theorem ortho_03 : krawtchouk_inner 0 3 = 0 := by native_decide
theorem ortho_04 : krawtchouk_inner 0 4 = 0 := by native_decide
theorem ortho_12 : krawtchouk_inner 1 2 = 0 := by native_decide
theorem ortho_13 : krawtchouk_inner 1 3 = 0 := by native_decide
theorem ortho_14 : krawtchouk_inner 1 4 = 0 := by native_decide
theorem ortho_23 : krawtchouk_inner 2 3 = 0 := by native_decide
theorem ortho_24 : krawtchouk_inner 2 4 = 0 := by native_decide
theorem ortho_34 : krawtchouk_inner 3 4 = 0 := by native_decide

/-- Summary: Krawtchouk orthogonality holds for all (w,w') pairs at n=4. -/
theorem krawtchouk_orthogonality_n4 :
    -- Off-diagonal vanishing
    krawtchouk_inner 0 1 = 0 ∧ krawtchouk_inner 0 2 = 0 ∧
    krawtchouk_inner 0 3 = 0 ∧ krawtchouk_inner 0 4 = 0 ∧
    krawtchouk_inner 1 2 = 0 ∧ krawtchouk_inner 1 3 = 0 ∧
    krawtchouk_inner 1 4 = 0 ∧ krawtchouk_inner 2 3 = 0 ∧
    krawtchouk_inner 2 4 = 0 ∧ krawtchouk_inner 3 4 = 0 ∧
    -- Diagonal = 16 · C(4,w)
    krawtchouk_inner 0 0 = 16 * 1 ∧ krawtchouk_inner 1 1 = 16 * 4 ∧
    krawtchouk_inner 2 2 = 16 * 6 ∧ krawtchouk_inner 3 3 = 16 * 4 ∧
    krawtchouk_inner 4 4 = 16 * 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> native_decide
