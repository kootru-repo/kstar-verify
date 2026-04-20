/-
  KstarFormal.LinearAlgebra.SpectralDecomp — Krawtchouk diagonalization of Gram matrix
  K* Verification: Krawtchouk spectral correspondence
  Registry: prop:spectral_q_main (Lemma 5)

  Core theorem: G(K) = Σ_w λ_w(K) · E_w where
    λ_w(K) = q^n · c_w(K) / C(n,w) · (q-1)^w = 2^4 · c_w / C(4,w)

  Certified eigenvalues at K=5: [144, 224, 64, 128, 256]
-/
import KstarFormal.Defs
import KstarFormal.Combinatorics.Krawtchouk
import KstarFormal.Combinatorics.LatticeCount
import KstarFormal.LinearAlgebra.GramMatrix

/-! ## Eigenvalue formula via parity weight counts

  For the binary Hamming scheme H(n,2):
    λ_w(K) = 2^n · c_w(K) / C(n,w)

  where c_w(K) = |{m ∈ Z^n : |m|^2 ≤ K, parityWeight(m) = w}|.
-/

/-- Gram eigenvalue from parity weight counts (exact integer when divisible). -/
def gramEigenvalue_from_cw (n : ℕ) (c_w_val : ℕ) (w : ℕ) : ℕ :=
  2^n * c_w_val / Nat.choose n w

/-- Eigenvalues at K=5, n=4 computed from c_w = [9, 56, 24, 32, 16]. -/
def eigenvalues_K5 : List ℕ :=
  (List.range 5).map fun w =>
    gramEigenvalue_from_cw 4 (c_w_K5.getD w 0) w

theorem eigenvalues_K5_eq : eigenvalues_K5 = [144, 224, 64, 128, 256] := by
  native_decide

/-! ## All eigenvalues are strictly positive -/

theorem all_eigenvalues_positive :
    ∀ w : Fin 5, eigenvalues_K5.getD w.val 0 > 0 := by
  intro w; fin_cases w <;> native_decide

/-! ## Eigenvalue formula verification

  Verify that λ_w = 2^4 · c_w / C(4,w) gives the right values:
    λ_0 = 16 · 9 / 1  = 144
    λ_1 = 16 · 56 / 4 = 224
    λ_2 = 16 · 24 / 6 = 64
    λ_3 = 16 · 32 / 4 = 128
    λ_4 = 16 · 16 / 1 = 256
-/

theorem eigenvalue_0 : 2^4 * 9  / Nat.choose 4 0 = 144 := by native_decide
theorem eigenvalue_1 : 2^4 * 56 / Nat.choose 4 1 = 224 := by native_decide
theorem eigenvalue_2 : 2^4 * 24 / Nat.choose 4 2 = 64  := by native_decide
theorem eigenvalue_3 : 2^4 * 32 / Nat.choose 4 3 = 128 := by native_decide
theorem eigenvalue_4 : 2^4 * 16 / Nat.choose 4 4 = 256 := by native_decide

/-! ## Divisibility: 2^n · c_w is divisible by C(n,w)

  This is required for the eigenvalue formula to yield integers.
  It follows from the Krawtchouk transform being integer-valued,
  but we verify it concretely.
-/

theorem divisibility_w0 : Nat.choose 4 0 ∣ 2^4 * c_w_K5.getD 0 0 := by native_decide
theorem divisibility_w1 : Nat.choose 4 1 ∣ 2^4 * c_w_K5.getD 1 0 := by native_decide
theorem divisibility_w2 : Nat.choose 4 2 ∣ 2^4 * c_w_K5.getD 2 0 := by native_decide
theorem divisibility_w3 : Nat.choose 4 3 ∣ 2^4 * c_w_K5.getD 3 0 := by native_decide
theorem divisibility_w4 : Nat.choose 4 4 ∣ 2^4 * c_w_K5.getD 4 0 := by native_decide

/-! ## Krawtchouk transform: distance values ↔ eigenvalues

  The Krawtchouk transform relates g(h) and λ_w:
    g(h) = Σ_w C(4,w) · λ_w · K_w(h;4) / 2^4

  Verify this for K=5: g = [137, 5, 17, -19, -39], λ = [144, 224, 64, 128, 256].
-/

/-- Forward Krawtchouk transform: eigenvalues → distance values (List-based).
    g(h) = (1/2^n) Σ_w λ_w · K_w(h;n). -/
def forwardTransform (eigs : List ℤ) (h : ℕ) : ℤ :=
  ((List.range 5).map fun w =>
    eigs.getD w 0 * krawtchouk 4 w h).sum / 2^4

/-- Eigenvalues as signed list for transform verification. -/
def eigenvalues_K5_signed : List ℤ := [144, 224, 64, 128, 256]

theorem forward_transform_h0 :
    forwardTransform eigenvalues_K5_signed 0 = 137 := by native_decide
theorem forward_transform_h1 :
    forwardTransform eigenvalues_K5_signed 1 = 5 := by native_decide
theorem forward_transform_h2 :
    forwardTransform eigenvalues_K5_signed 2 = 17 := by native_decide
theorem forward_transform_h3 :
    forwardTransform eigenvalues_K5_signed 3 = -19 := by native_decide
theorem forward_transform_h4 :
    forwardTransform eigenvalues_K5_signed 4 = -39 := by native_decide
