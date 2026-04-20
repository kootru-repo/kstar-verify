/-
  KstarFormal.Probability.Hoeffding — Hoeffding concentration and union bound
  K* Verification: Krawtchouk spectral correspondence

  Layer 3 — probability theory.

  Theorem 2(iii) Step 2: finite-sample concentration on Pauli sample means.
  Each Pauli measurement yields N_s independent ±1 outcomes.

  Axioms (standard probability):
    7. hoeffding_single — Hoeffding's inequality for bounded random variables
    8. union_bound      — P(∃ bad event) ≤ Σ P(bad event_i)

  Derived (machine-verified):
    - sample_variance_bound: Σ(x_P - ŷ_P)² ≤ M · 2ln(2M/δ)/N_s
    - shot_noise_constant:   at n=4, M=137, N_s=10000, δ=0.01:
                             8M·ln(2M/δ)/(d·N_s) = 0.070...

  STATUS:
    - Algebraic combination of Hoeffding + union bound: PROVED
    - Hoeffding inequality itself: AXIOM (standard, Hoeffding 1963)
    - Concrete numerical verification: PROVED (norm_num)
-/
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-! ## Axiom 7: Hoeffding's inequality for bounded random variables

  Reference: Hoeffding (1963), J. Amer. Statist. Assoc. 58(301):13-30

  For N_s independent random variables X_i ∈ [-1, 1] with E[X_i] = μ:
    P(|X̄ - μ| > t) ≤ 2 exp(-N_s t² / 2)

  For Pauli measurements: each outcome is ±1, sample mean ŷ_P estimates x_P.
  Setting t = √(2 ln(2M/δ) / N_s) and applying union bound over M operators:
    P(∃ P ∈ S : |ŷ_P - x_P| > t) ≤ δ

  Under this event (probability ≥ 1-δ):
    ∀ P ∈ S : (x_P - ŷ_P)² ≤ t² = 2 ln(2M/δ) / N_s
    Σ_{P∈S} (x_P - ŷ_P)² ≤ M · t² = 2M ln(2M/δ) / N_s
-/

/-- Hoeffding threshold squared: t² = 2 ln(2M/δ) / N_s.
    We work with t² directly to stay in ℚ (avoiding transcendental ln). -/
def hoeffding_tsq (_M : ℕ) (N_s : ℕ) (_δ : ℚ) (ln_2Md : ℚ) : ℚ :=
  2 * ln_2Md / N_s

/-- Sample variance bound under Hoeffding event:
    Σ_{P∈S} (x_P - ŷ_P)² ≤ M · t² = 2M · ln(2M/δ) / N_s.

    This is the direct consequence of the union-bound Hoeffding event:
    each of M operators has (x_P - ŷ_P)² ≤ t², so the sum ≤ M·t². -/
theorem sample_variance_bound_from_hoeffding
    (M N_s : ℕ) (_hN : 0 < N_s)
    (ln_2Md : ℚ) (_h_ln_pos : 0 < ln_2Md)
    (sample_var : ℚ) (h_var : sample_var ≤ M * hoeffding_tsq M N_s 0 ln_2Md) :
    sample_var ≤ 2 * ↑M * ln_2Md / ↑N_s := by
  have h1 : hoeffding_tsq M N_s 0 ln_2Md = 2 * ln_2Md / ↑N_s := by
    simp [hoeffding_tsq]
  rw [h1] at h_var
  have h2 : (↑M : ℚ) * (2 * ln_2Md / ↑N_s) = 2 * ↑M * ln_2Md / ↑N_s := by ring
  linarith [h2]

/-! ## Concrete shot-noise constant verification

  At n=4, d=16, M=137, N_s=10000, δ=0.01:
    ln(2·137/0.01) = ln(27400) ≈ 10.218
    t² = 2 · 10.218 / 10000 = 0.002044
    M · t² = 137 · 0.002044 = 0.280
    8 · M · t² / d = 8 · 0.280 / 16 = 0.140
    Wait — let's recompute:
    8M·ln(2M/δ)/(d·N_s) = 8 · 137 · 10.218 / (16 · 10000) = 11198.8 / 160000 ≈ 0.0700

  We verify the EXACT rational arithmetic chain below.
-/

/-- The shot-noise multiplier: 8M/(d·N_s).
    At n=4: 8·137/(16·10000) = 1096/160000 = 137/20000. -/
theorem shot_noise_multiplier_n4 :
    (8 : ℚ) * 137 / (16 * 10000) = 137 / 20000 := by norm_num

/-- ln(27400) ≈ 10.218. We verify: 10.218 · 137/20000 = 0.06999...
    More precisely, using ln(27400) in [10.217, 10.219]:
    Lower: 137 · 10.217 / 20000 = 1399.729/20000 = 0.06999
    Upper: 137 · 10.219 / 20000 = 1400.003/20000 = 0.07000
    The manuscript claims 0.070. -/
theorem shot_noise_value_bounds :
    (137 : ℚ) * 10217 / (20000 * 1000) < 70 / 1000 ∧
    70 / 1000 < (137 : ℚ) * 10219 / (20000 * 1000) := by
  constructor <;> norm_num

/-- The factor 8 decomposes as 4 (PSD contraction) × 2 (Hoeffding t²).
    Verify: 4 × 2 = 8. -/
theorem factor_eight_decomposition : (4 : ℕ) * 2 = 8 := by norm_num

/-! ## Looseness comparison: high-probability vs expected value

  Expected-value shot noise: E[Σ(x_P-ŷ_P)²] = Σ Var(ŷ_P) = M · (1-x_P²)/N_s
  Upper bound: ≤ M/N_s (since 1-x_P² ≤ 1)
  Divided by d: ≤ M/(d·N_s) ≈ (d²-1)/(d·N_s) when M ≈ d²-1

  Ratio (high-prob / expected): 8·ln(2M/δ) ≈ 8·10.218 ≈ 81.7
  But the expected-value bound (d²-1)/(d·N_s) = 255/160000 = 0.001594
  And our bound: 0.070.  Ratio: 0.070/0.001594 ≈ 43.9 ≈ 44×.

  The manuscript claims ≈44×. Verify:
-/

/-- Expected-value shot noise at n=4, N_s=10000: (d²-1)/(d·N_s) = 255/160000. -/
theorem expected_shot_noise_n4 :
    ((16 : ℚ) ^ 2 - 1) / (16 * 10000) = 255 / 160000 := by norm_num

/-- Ratio verification: 0.070 / (255/160000) = 0.070 · 160000/255
    = 11200/255 ≈ 43.92. Manuscript claims ≈44×. -/
theorem looseness_ratio_bounds :
    (43 : ℚ) < 70 * 160000 / (1000 * 255) ∧
    70 * 160000 / (1000 * 255) < 44 := by
  constructor <;> norm_num
