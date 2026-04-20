/-
  KstarFormal.Certified — Layer 1 complete assembly
  K* Verification: Krawtchouk spectral correspondence

  This file imports all Layer 1 (pure combinatorics) results and re-exports
  the key theorems. If this file compiles with no axiom gaps, the combinatorial
  core of the paper is machine-verified.

  Machine-verified claims (Layer 1):
    ✓ N_4(5) = 137 (lattice count)
    ✓ c_w(K=5) = [9, 56, 24, 32, 16] (parity weight distribution)
    ✓ λ_w(K=5) = [144, 224, 64, 128, 256] (Gram eigenvalues)
    ✓ All λ_w > 0 (Gram matrix positive definite)
    ✓ κ(G) = 256/64 = 4 (condition number)
    ✓ M_w(K=5) = [1, 12, 54, 54, 16] (greedy redistribution)
    ✓ w_sat = 2 (weight saturation index)
    ✓ K* = 5 (minimal saturation cutoff)
    ✓ |S| ≥ 66 (operator lower bound, Corollary 2)
    ✓ λ_w(K=4) ≤ λ_w(K=5) (monotonicity, Lemma 3)
    ✓ Krawtchouk orthogonality verified for n=4

  Formal statements with axiom gaps (Layer 2/3, not yet verified):
    ○ Lemma 1 (Hessian diagonal structure)
    ○ Lemma 2 (purity bound on positivity excursion)
    ○ Theorem 1 (basin separation)
    ○ Theorem 2 (spectral characterization — general n)
    ○ Corollary 1 (approximate locality)
    ○ Lemma 4 (hypergeometric bound)
-/

-- Layer 1: Pure combinatorics (gap-free target)
import KstarFormal.Combinatorics.WeightSaturation
import KstarFormal.LinearAlgebra.Eigenvalues
import KstarFormal.LinearAlgebra.Monotonicity

/-!
## Layer 1 Certificate

All theorems below are re-exports from the imported files.
They are listed here as the formal certificate document.
-/

-- 1. Lattice count
#check @N4_five           -- N_4(5) = 137
#check @c_w_K5_sum        -- sum(c_w) = 137

-- 2. Eigenvalue computation
#check @eigenvalues_K5_eq -- λ = [144, 224, 64, 128, 256]
#check @all_eigenvalues_positive  -- ∀ w, λ_w > 0

-- 3. Greedy redistribution
#check @M_w_K5_eq         -- M_w = [1, 12, 54, 54, 16]
#check @M_w_K5_sum        -- sum(M_w) = 137

-- 4. Weight saturation
#check @kstar_eq_five     -- K* = 5
#check @w_sat_K5_correct  -- w_sat = 2

-- 5. Lower bound (Corollary 2)
#check @operator_lower_bound_eq  -- bound = 66
#check @kstar_exceeds_lower_bound  -- 137 > 66

-- 6. Monotonicity (Lemma 3)
#check @eigenvalue_monotone_K4_K5  -- λ_w(4) ≤ λ_w(5)
#check @M_w_monotone_K4_K5        -- M_w(4) ≤ M_w(5)

-- 7. Support-completeness (Lemma 6)
#check @support_completeness_n4    -- all λ_w > 0 → all c_w > 0
#check @full_rank_K5               -- rank G(5) = 16

-- 8. Krawtchouk orthogonality
#check @krawtchouk_orthogonality_n4

-- 9. Condition number
#check @condition_number            -- κ = 4
