/-
  KstarFormal.Quantum.HradilDerivation — Hradil R-operator closed form (Tier 1C)
  K* Verification: Krawtchouk spectral correspondence

  Formalizes the algebraic identity behind the corrected Hradil R-operator
  for the Pauli-binomial likelihood used in robust MLE state reconstruction.

  For a single Pauli measurement i with outcome bit b (= ±1 frequency
  encoded as a rational in [-1, 1]) and current model expectation y
  (also in (-1, 1)), the per-measurement Hradil contribution sums the
  POVM-weighted update over the two outcomes:

      R_pre(b, y) = (1+b)/(1+y) · (I + P)/4
                  + (1-b)/(1-y) · (I - P)/4

  where I and P are formal generators of a free 2-dimensional ℚ-module.
  Direct algebraic simplification (clearing the common denominator
  (1 - y²) and combining the (I, P) coefficients) gives the pole-free
  closed form

      R_closed(b, y) = (1 - b·y) / (1 - y²) · (I/2)
                     + (b - y) / (1 - y²) · (P/2).

  This eliminates the (1 ± y) poles that the b/y approximation form
  carries, removing the need for ratio clipping in the iteration.

  The proof works coordinate-wise in the free module ℚ × ℚ
  (with I = (1, 0) and P = (0, 1)), reducing the equality of the two
  R-operators to two scalar identities over ℚ that are dispatched by
  `field_simp` + `ring`.

  Status: Tier 1C formalization, sorry-free. Treats I and P as formal
  generators (no Hilbert-space structure required); the manuscript-level
  identity is the algebraic content.
-/
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

namespace KstarFormal.HradilDerivation

/-! ## Per-measurement Hradil contribution: pre- and closed forms

  We work in the free ℚ-module ℚ × ℚ with the convention
    I := (1, 0)    "identity component"
    P := (0, 1)    "Pauli component"
  so an element (a, c) of ℚ × ℚ stands for `a · I + c · P`.
-/

/-- Pre-simplification Hradil contribution as a (coeffI, coeffP) pair.
    Encodes
      ((1+b)/(1+y)) · (I + P)/4  +  ((1-b)/(1-y)) · (I - P)/4
    componentwise in the free ℚ-module ℚ × ℚ. -/
def R_pre (b y : ℚ) : ℚ × ℚ :=
  ( (1 + b) / (1 + y) / 4 + (1 - b) / (1 - y) / 4
  , (1 + b) / (1 + y) / 4 - (1 - b) / (1 - y) / 4 )

/-- Closed-form Hradil contribution after clearing the common denominator.
    Encodes
      (1 - b·y)/(1 - y²) · (I/2)  +  (b - y)/(1 - y²) · (P/2). -/
def R_closed (b y : ℚ) : ℚ × ℚ :=
  ( (1 - b * y) / (1 - y ^ 2) / 2
  , (b - y)     / (1 - y ^ 2) / 2 )

/-! ## The closed-form identity -/

/-- **Tier 1C main lemma — coefficient of I.** The identity-component of
    the per-measurement Hradil contribution simplifies to (1 - by) / (2(1 - y²)). -/
theorem R_pre_I_eq (b y : ℚ) (hpos : 1 + y ≠ 0) (hneg : 1 - y ≠ 0) :
    (R_pre b y).1 = (R_closed b y).1 := by
  unfold R_pre R_closed
  have hsq : 1 - y ^ 2 ≠ 0 := by
    have : (1 - y ^ 2) = (1 + y) * (1 - y) := by ring
    rw [this]
    exact mul_ne_zero hpos hneg
  field_simp
  ring

/-- **Tier 1C main lemma — coefficient of P.** The Pauli-component of
    the per-measurement Hradil contribution simplifies to (b - y) / (2(1 - y²)). -/
theorem R_pre_P_eq (b y : ℚ) (hpos : 1 + y ≠ 0) (hneg : 1 - y ≠ 0) :
    (R_pre b y).2 = (R_closed b y).2 := by
  unfold R_pre R_closed
  have hsq : 1 - y ^ 2 ≠ 0 := by
    have : (1 - y ^ 2) = (1 + y) * (1 - y) := by ring
    rw [this]
    exact mul_ne_zero hpos hneg
  field_simp
  ring

/-- **Tier 1C main theorem.** The full per-measurement R-operator
    contribution agrees in both forms whenever y ≠ ±1 (the only
    constraint of the Pauli expectation regime, automatic for any
    quantum state strictly inside the Bloch ball). -/
theorem R_pre_eq_R_closed (b y : ℚ) (hpos : 1 + y ≠ 0) (hneg : 1 - y ≠ 0) :
    R_pre b y = R_closed b y := by
  ext
  · exact R_pre_I_eq b y hpos hneg
  · exact R_pre_P_eq b y hpos hneg

/-! ## Pole-freeness at y = 0

  The closed form is well-defined at y = 0 (where the pre-form is
  vacuously well-defined, but the b/y approximation that the manuscript
  replaces would already exhibit a 0/0 issue). At y = 0, the closed
  form reduces to (I + b · P) / 2, which is the maximally mixed state
  perturbed in the b direction — exactly what the Hradil iteration
  should produce on the first step from R = I/d.
-/

/-- At y = 0, the closed form reduces to (I + b·P)/2. -/
theorem R_closed_at_zero (b : ℚ) :
    R_closed b 0 = (1 / 2, b / 2) := by
  unfold R_closed
  refine Prod.mk.injEq _ _ _ _ |>.mpr ⟨?_, ?_⟩ <;> norm_num

/-! ## Fixed-point behaviour

  At b = y (i.e. the model already matches the data), the closed
  form should give (1/(2(1-y²)) · (1 - y²) · I/1, 0) = I/2, the
  identity (modulo the per-measurement normalization 1/2) — confirming
  that the iteration is at a fixed point.
-/

/-- At b = y the Pauli component vanishes. -/
theorem R_closed_fixed_point_P (y : ℚ) (hpos : 1 + y ≠ 0) (hneg : 1 - y ≠ 0) :
    (R_closed y y).2 = 0 := by
  unfold R_closed
  have hsq : 1 - y ^ 2 ≠ 0 := by
    have : (1 - y ^ 2) = (1 + y) * (1 - y) := by ring
    rw [this]; exact mul_ne_zero hpos hneg
  rw [show (y - y : ℚ) = 0 from by ring]
  simp

/-- At b = y the identity component is exactly 1/2 (the per-measurement
    contribution to the identity, summing over two outcomes). -/
theorem R_closed_fixed_point_I (y : ℚ) (hpos : 1 + y ≠ 0) (hneg : 1 - y ≠ 0) :
    (R_closed y y).1 = 1 / 2 := by
  unfold R_closed
  have hsq : 1 - y ^ 2 ≠ 0 := by
    have : (1 - y ^ 2) = (1 + y) * (1 - y) := by ring
    rw [this]; exact mul_ne_zero hpos hneg
  field_simp

/-- **Bundled fixed-point characterization.** When the model expectation
    matches the data (`b = y`), the per-measurement contribution to the
    Hradil R-operator collapses to the pure-identity element `(1/2, 0)`
    in the (I, P) basis: identity coefficient 1/2, Pauli coefficient 0.
    The R operator is therefore proportional to the identity at the
    fixed point, so the iteration ρ ↦ R ρ R / tr(R ρ R) is stationary. -/
theorem R_closed_fixed_point (y : ℚ) (hpos : 1 + y ≠ 0) (hneg : 1 - y ≠ 0) :
    R_closed y y = (1 / 2, 0) := by
  ext
  · exact R_closed_fixed_point_I y hpos hneg
  · exact R_closed_fixed_point_P y hpos hneg

/-! ## Concrete numerical anchors and negative witnesses

  Spot-check the closed form on hard-coded (b, y) pairs cross-validated
  against the Python harness (`cross_validate_lean.py`, Section 7), plus
  a negative anchor proving the closed form is *not* identically zero
  (rules out the trivial-degeneration class of interpretation errors).
-/

/-- Numerical anchor: at (b, y) = (1/2, 1/3) the closed form gives
    coefficients (15/32, 3/32). Cross-checked against Python. -/
theorem R_closed_anchor_half_third :
    R_closed (1/2) (1/3) = (15/32, 3/32) := by
  unfold R_closed; refine Prod.mk.injEq _ _ _ _ |>.mpr ⟨?_, ?_⟩ <;> norm_num

/-- Numerical anchor at (b, y) = (-1/4, 1/2). -/
theorem R_closed_anchor_neg :
    R_closed (-1/4) (1/2) = (3/4, -1/2) := by
  unfold R_closed; refine Prod.mk.injEq _ _ _ _ |>.mpr ⟨?_, ?_⟩ <;> norm_num

/-- Negative anchor: when b ≠ y the Pauli coefficient is nonzero —
    the iteration moves off the fixed point in the right direction. -/
theorem R_closed_moves_when_b_ne_y :
    (R_closed (1/2) (0 : ℚ)).2 ≠ 0 := by
  unfold R_closed; norm_num

/-- Negative anchor: the closed form is *not* the trivial map sending
    everything to (0, 0). At (b, y) = (1, 0), the I-coefficient is 1/2. -/
theorem R_closed_nontrivial :
    (R_closed 1 (0 : ℚ)).1 ≠ 0 := by
  unfold R_closed; norm_num

end KstarFormal.HradilDerivation
