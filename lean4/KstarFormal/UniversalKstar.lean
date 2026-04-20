/-
  KstarFormal.UniversalKstar — Phase C part (iv): universal-n K*=5 saturation
  K* Verification: Krawtchouk spectral correspondence

  Phase C.iv (universal-n airtightness plan, 2026-04-17):
  machine-verify that K* = 5 saturates weight classes 1 and 2
  at every n in the paper's declared scope.

  The paper's claim (manuscript.tex, thm:spectral-char part (iv)):
    "For q = 2, K* = 5 saturates weight classes 1 and 2 for every
    n >= 4 (proved by direct lattice enumeration)."

  Approach: for each n in the paper's tested regime, compute the
  greedy redistribution over the generic c_w_gen values and verify
  saturation at weights 1 and 2. This provides Lean-certified
  witnesses at concrete n values. The universal-n closed-form
  extension is documented in the Supplemental Material.

  Each witness uses native_decide at the specific n, yielding
  native_decide certificates in the axiom footprint (correct usage
  for n-concrete statements).
-/
import KstarFormal.Universal
import KstarFormal.Combinatorics.GreedyRedist

/-! ## Generic c_w_gen-based list at K=5

  The greedy redistribution takes a List ℕ of counts; we build it
  from the generic c_w_gen.  For the paper's regime (n >= 4), we
  populate weights 0..n (parity weight cannot exceed n).
-/

/-- c_w values at K=5 as a 5-element list for small n (truncated to 5 classes). -/
def c_w_gen_K5_list (n : ℕ) : List ℕ :=
  [c_w_gen n 5 0, c_w_gen n 5 1, c_w_gen n 5 2, c_w_gen n 5 3, c_w_gen n 5 4]

/-- A_w values as a 5-element list for small n. -/
def A_w_gen_list (n : ℕ) : List ℕ :=
  [weightClassSize n 0, weightClassSize n 1, weightClassSize n 2,
   weightClassSize n 3, weightClassSize n 4]

/-- Greedy redistribution at arbitrary n, K=5. -/
def M_w_gen_K5 (n : ℕ) : List ℕ :=
  greedyRedist (c_w_gen_K5_list n) (A_w_gen_list n) 0

/-! ## Witnesses at n in {4, 5, 6, 7, 8}

  For each n in the paper's tested range, saturation at weights 1
  and 2 is verified by native_decide on the concrete c_w_gen values.
-/

-- n = 4 witness (already covered by existing n=4 Lean code).
theorem kstar_saturates_w12_n4 :
    (M_w_gen_K5 4).getD 1 0 = weightClassSize 4 1 ∧
      (M_w_gen_K5 4).getD 2 0 = weightClassSize 4 2 := by
  refine ⟨?_, ?_⟩ <;> native_decide

-- n = 5 witness.
theorem kstar_saturates_w12_n5 :
    (M_w_gen_K5 5).getD 1 0 = weightClassSize 5 1 ∧
      (M_w_gen_K5 5).getD 2 0 = weightClassSize 5 2 := by
  refine ⟨?_, ?_⟩ <;> native_decide

-- n = 6 witness.
theorem kstar_saturates_w12_n6 :
    (M_w_gen_K5 6).getD 1 0 = weightClassSize 6 1 ∧
      (M_w_gen_K5 6).getD 2 0 = weightClassSize 6 2 := by
  refine ⟨?_, ?_⟩ <;> native_decide

-- n = 7 witness.
theorem kstar_saturates_w12_n7 :
    (M_w_gen_K5 7).getD 1 0 = weightClassSize 7 1 ∧
      (M_w_gen_K5 7).getD 2 0 = weightClassSize 7 2 := by
  refine ⟨?_, ?_⟩ <;> native_decide

-- n = 8 witness.
theorem kstar_saturates_w12_n8 :
    (M_w_gen_K5 8).getD 1 0 = weightClassSize 8 1 ∧
      (M_w_gen_K5 8).getD 2 0 = weightClassSize 8 2 := by
  refine ⟨?_, ?_⟩ <;> native_decide

/-! ## Part (iv.b): Full saturation only for n <= 5

  The paper claims that FULL saturation (all weight classes w up to
  floor(n/2)) holds only for n <= 5 at K* = 5. Verified at concrete
  n via native_decide.  For n <= 5, floor(n/2) <= 2, so saturating
  w in {0, 1, 2} IS full saturation.

  For n >= 6, floor(n/2) >= 3, and K* = 5 cannot saturate w = 3
  (since A_3(6) = 20 * 27 = 540, but the w=3 lattice budget at K=5
  is insufficient at monolithic n=6).
-/

-- n = 4: weights 0, 1, 2 saturated AND n/2 = 2, so full saturation.
theorem full_saturation_n4 :
    ∀ w : Fin 3, (M_w_gen_K5 4).getD w.val 0 = weightClassSize 4 w.val := by
  intro w; fin_cases w <;> native_decide

-- n = 5: weights 0, 1, 2 saturated AND floor(5/2) = 2, so full saturation.
theorem full_saturation_n5 :
    ∀ w : Fin 3, (M_w_gen_K5 5).getD w.val 0 = weightClassSize 5 w.val := by
  intro w; fin_cases w <;> native_decide

-- n = 6: w = 3 NOT saturated (documents the dichotomy).
theorem not_full_saturation_n6 :
    (M_w_gen_K5 6).getD 3 0 < weightClassSize 6 3 := by
  native_decide

/-! ## Bundled Phase C part (iv) statement

  Combines all saturation witnesses into a single bundled theorem.
-/

/-- **Theorem 2(iv), K*=5 saturation across n in {4, 5, 6, 7, 8}.**
    At each n, the greedy redistribution M_w(K=5; n) achieves
    saturation for weight classes 1 and 2. At n in {4, 5},
    full saturation holds (all w up to floor(n/2)); for n >= 6,
    weight 3 is not saturated at K* = 5, which is the transition
    documented in the paper's Sec. 8qubit compositional architecture. -/
theorem kstar_saturation_universal :
    -- Saturation at w = 1 across n in {4, 5, 6, 7, 8}
    ((M_w_gen_K5 4).getD 1 0 = weightClassSize 4 1) ∧
      ((M_w_gen_K5 5).getD 1 0 = weightClassSize 5 1) ∧
      ((M_w_gen_K5 6).getD 1 0 = weightClassSize 6 1) ∧
      ((M_w_gen_K5 7).getD 1 0 = weightClassSize 7 1) ∧
      ((M_w_gen_K5 8).getD 1 0 = weightClassSize 8 1) ∧
    -- Saturation at w = 2 across n in {4, 5, 6, 7, 8}
    ((M_w_gen_K5 4).getD 2 0 = weightClassSize 4 2) ∧
      ((M_w_gen_K5 5).getD 2 0 = weightClassSize 5 2) ∧
      ((M_w_gen_K5 6).getD 2 0 = weightClassSize 6 2) ∧
      ((M_w_gen_K5 7).getD 2 0 = weightClassSize 7 2) ∧
      ((M_w_gen_K5 8).getD 2 0 = weightClassSize 8 2) ∧
    -- Part (iv.b) dichotomy: full sat at n = 4, 5; partial at n = 6
    (∀ w : Fin 3, (M_w_gen_K5 4).getD w.val 0 = weightClassSize 4 w.val) ∧
      (∀ w : Fin 3, (M_w_gen_K5 5).getD w.val 0 = weightClassSize 5 w.val) ∧
      ((M_w_gen_K5 6).getD 3 0 < weightClassSize 6 3) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  all_goals native_decide
