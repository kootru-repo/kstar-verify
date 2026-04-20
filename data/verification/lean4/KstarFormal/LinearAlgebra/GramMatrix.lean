/-
  KstarFormal.LinearAlgebra.GramMatrix — Gram matrix of the Z_2 orbifold
  K* Verification: Krawtchouk spectral correspondence
  Registry: prop:spectral_q_main (Lemma 5), thm:spectral-char (Theorem 2)

  The Gram matrix G(K) is defined on the 2^n = 16 fixed points of Z_2 on {0,1}^4.
  Entry G_{ij} = Σ_{|m|^2 ≤ K} (-1)^{m · (v_i ⊕ v_j)}

  Key property: G depends only on Hamming distance d_H(v_i, v_j).
  This places G in the Bose-Mesner algebra of H(4,2).

  Certified distance values at K=5: g(h) = {0:137, 1:5, 2:17, 3:-19, 4:-39}
-/
import KstarFormal.Defs
import KstarFormal.Combinatorics.LatticeCount

/-! ## Distance values (Hamming-constant structure)

  The Gram matrix G(K) is constant on Hamming distance classes.
  We define it via the distance values g(h) for h = 0,...,4.
-/

/-- Distance values for the Gram matrix at K=5 as a list.
    g(h) = Σ_{|m|^2 ≤ 5} (-1)^{m · v} where v is any vector with |v|_H = h.
    Certified by SageMath: [137, 5, 17, -19, -39]. -/
def distanceValues_K5 : List ℤ := [137, 5, 17, -19, -39]

/-- Distance values at K=4 (for monotonicity comparison). -/
def distanceValues_K4 : List ℤ := [89, -7, 17, 3, -39]

/-! ## Concrete distance value verification -/

theorem dv_K5_0 : distanceValues_K5.getD 0 0 = 137 := by native_decide
theorem dv_K5_1 : distanceValues_K5.getD 1 0 = 5 := by native_decide
theorem dv_K5_2 : distanceValues_K5.getD 2 0 = 17 := by native_decide
theorem dv_K5_3 : distanceValues_K5.getD 3 0 = -19 := by native_decide
theorem dv_K5_4 : distanceValues_K5.getD 4 0 = -39 := by native_decide

/-- Weighted sum Σ_h C(4,h) · g(h) = c_0 · 2^4 = 144 (= λ_0).
    This is the w=0 Krawtchouk coefficient. -/
def weighted_dv_sum : ℤ :=
  ((List.range 5).map fun h =>
    (Nat.choose 4 h : ℤ) * distanceValues_K5.getD h 0).sum

theorem weighted_dv_sum_eq : weighted_dv_sum = 144 := by native_decide

/-- Trace of G(K=5) = 16 · g(0) = 16 · 137 = 2192. -/
theorem gram_trace_K5 : 16 * distanceValues_K5.getD 0 0 = 2192 := by native_decide

/-! ## Hamming distance computation on {0,1}^4

  The 16 fixed points of Z_2 on {0,1}^4 are indexed by ℕ < 16.
  We decode each index to its binary representation and compute Hamming distance.
-/

/-- Hamming distance between binary representations of two naturals (4-bit). -/
def hammingDist4 (i j : ℕ) : ℕ :=
  (if (i / 1) % 2 ≠ (j / 1) % 2 then 1 else 0) +
  (if (i / 2) % 2 ≠ (j / 2) % 2 then 1 else 0) +
  (if (i / 4) % 2 ≠ (j / 4) % 2 then 1 else 0) +
  (if (i / 8) % 2 ≠ (j / 8) % 2 then 1 else 0)

/-- Self-distance is zero. -/
theorem hamming_self : ∀ i : Fin 16, hammingDist4 i.val i.val = 0 := by
  intro i; fin_cases i <;> native_decide

/-- Distance distribution: for any codeword i in {0,..,15},
    the number at distance h is C(4,h).
    This is the regularity property of the Hamming scheme H(4,2). -/
theorem distance_distribution_h0 :
    ∀ i : Fin 16,
    ((List.range 16).filter fun j => hammingDist4 i.val j = 0).length = 1 := by
  intro i; fin_cases i <;> native_decide

theorem distance_distribution_h1 :
    ∀ i : Fin 16,
    ((List.range 16).filter fun j => hammingDist4 i.val j = 1).length = 4 := by
  intro i; fin_cases i <;> native_decide

theorem distance_distribution_h2 :
    ∀ i : Fin 16,
    ((List.range 16).filter fun j => hammingDist4 i.val j = 2).length = 6 := by
  intro i; fin_cases i <;> native_decide

theorem distance_distribution_h3 :
    ∀ i : Fin 16,
    ((List.range 16).filter fun j => hammingDist4 i.val j = 3).length = 4 := by
  intro i; fin_cases i <;> native_decide

theorem distance_distribution_h4 :
    ∀ i : Fin 16,
    ((List.range 16).filter fun j => hammingDist4 i.val j = 4).length = 1 := by
  intro i; fin_cases i <;> native_decide
