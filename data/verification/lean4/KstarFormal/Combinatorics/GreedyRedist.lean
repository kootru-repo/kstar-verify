/-
  KstarFormal.Combinatorics.GreedyRedist — Greedy redistribution algorithm
  K* Verification: Krawtchouk spectral correspondence
  Registry: lem:monotone (Lemma 3), thm:spectral-char (Theorem 2(iv))

  The greedy redistribution takes lattice counts c_w and weight class caps A_w,
  and produces the measurement allocation M_w by greedily filling each weight
  class, carrying excess forward.

  Key certified values:
    c_w(K=5) = [9, 56, 24, 32, 16]
    A_w(n=4) = [1, 12, 54, 108, 81]
    M_w(K=5) = [1, 12, 54, 54, 16]
-/
import KstarFormal.Defs
import KstarFormal.Combinatorics.LatticeCount

/-! ## Greedy redistribution algorithm -/

/-- Greedy redistribution: given counts `c` and caps `a`, produce allocation `m`.
    Processes weights in order, carrying forward any excess. -/
def greedyRedist : List ℕ → List ℕ → ℕ → List ℕ
  | [], _, _ => []
  | _ :: cs, [], e => 0 :: greedyRedist cs [] e  -- no cap remaining
  | c :: cs, a :: as, e =>
    let available := c + e
    let m := min available a
    let excess := available - m
    m :: greedyRedist cs as excess

/-! ## Concrete computation for n=4, K=5 -/

/-- M_w at K=5: greedy redistribution of c_w through A_w. -/
def M_w_K5 : List ℕ := greedyRedist c_w_K5 A_w_n4 0

theorem M_w_K5_eq : M_w_K5 = [1, 12, 54, 54, 16] := by native_decide

/-- Conservation: total allocation equals total count (no lattice points lost). -/
theorem M_w_K5_sum : M_w_K5.sum = 137 := by native_decide

/-- Weight saturation: weights 0, 1, 2 are fully saturated (M_w = A_w). -/
theorem weight_0_saturated : M_w_K5.getD 0 0 = A_w_n4.getD 0 0 := by native_decide
theorem weight_1_saturated : M_w_K5.getD 1 0 = A_w_n4.getD 1 0 := by native_decide
theorem weight_2_saturated : M_w_K5.getD 2 0 = A_w_n4.getD 2 0 := by native_decide

/-- Weight 3 is NOT saturated: M_3 = 54 < A_3 = 108. -/
theorem weight_3_not_saturated : M_w_K5.getD 3 0 < A_w_n4.getD 3 0 := by native_decide

/-- The saturation index w_sat = 2 (highest w where M_w = A_w). -/
def w_sat_K5 : ℕ := 2

theorem w_sat_K5_correct :
    (M_w_K5.getD w_sat_K5 0 = A_w_n4.getD w_sat_K5 0) ∧
    (M_w_K5.getD (w_sat_K5 + 1) 0 < A_w_n4.getD (w_sat_K5 + 1) 0) := by
  constructor <;> native_decide

/-! ## M_w at K=4 (for monotonicity comparison) -/

def c_w_K4 : List ℕ := [9, 8, 24, 32, 16]

theorem c_w_K4_sum : c_w_K4.sum = 89 := by native_decide

def M_w_K4 : List ℕ := greedyRedist c_w_K4 A_w_n4 0

theorem M_w_K4_eq : M_w_K4 = [1, 12, 28, 32, 16] := by native_decide

/-! ## Monotonicity: M_w(K=4) ≤ M_w(K=5) pointwise.
    This is a concrete instance of Lemma 3. -/

theorem M_w_monotone_K4_K5 :
    ∀ w : Fin 5, M_w_K4.getD w.val 0 ≤ M_w_K5.getD w.val 0 := by
  intro w; fin_cases w <;> native_decide
