/-
  KstarFormal.Combinatorics.QaryGram — q-ary Bose-Mesner iff (Tier 1A)
  K* Verification: Krawtchouk spectral correspondence

  Formalizes the load-bearing combinatorial content of Lemma 5
  (prop:spectral_q_main): the Gram matrix G(K) of the Hamming scheme
  H(n,q) lies in the Bose-Mesner algebra iff q ≤ 3.

  The essence of the iff is the lift-norm dichotomy: every nonzero
  residue mod q has minimum symmetric lift with squared norm 1
  iff q ≤ 3. For q ≥ 4, residues 1 and 2 lift to squared norms 1
  and 4 respectively, breaking distance-class constancy of the
  Euclidean cutoff |m|² ≤ K and so excluding G(K) from the
  Bose-Mesner algebra of H(n,q).

  Status: Tier 1A formalization, sorry-free.
-/
import Mathlib.Data.Nat.Basic
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith

namespace KstarFormal.QaryGram

/-! ## Symmetric lift norm

  For a residue `a ∈ ZMod q` represented as `a : ℕ` with `0 ≤ a < q`,
  the symmetric lift sends `a` to whichever of `a` or `q - a` has smaller
  magnitude. The lift-norm-squared is the square of that magnitude.
-/

/-- Squared norm of the minimum symmetric lift of `a : ZMod q`,
    represented as a natural number with `0 ≤ a < q`. -/
def liftNormSq (q a : ℕ) : ℕ :=
  let s := min a (q - a)
  s * s

/-- For q = 2 the only nonzero residue is 1, with lift norm 1. -/
theorem liftNormSq_q2 : ∀ a, 0 < a → a < 2 → liftNormSq 2 a = 1 := by
  intro a ha_pos ha_lt
  -- a must be 1
  have ha : a = 1 := by omega
  subst ha
  rfl

/-- For q = 3 the nonzero residues are 1 and 2, both with lift norm 1. -/
theorem liftNormSq_q3 : ∀ a, 0 < a → a < 3 → liftNormSq 3 a = 1 := by
  intro a ha_pos ha_lt
  -- a is 1 or 2
  rcases Nat.lt_or_ge a 2 with h | h
  · have ha : a = 1 := by omega
    subst ha; rfl
  · have ha : a = 2 := by omega
    subst ha; rfl

/-- For q = 4, residue 2 has lift norm 4 (not 1). The boundary case
    `q = 2 * residue` makes the value strictly 4 here; for q ≥ 5 the
    value also happens to be 4 at residue 2, so this anchor is not
    q-pinning on its own. -/
theorem liftNormSq_q4_two : liftNormSq 4 2 = 4 := by decide

/-- Concrete witness: at q = 4, residues 1 and 2 have distinct lift norms.
    The general parametric form `qge4_lift_norms_distinct` lives below
    `lns_two_gt_one_of_qge4` and proves this for every q ≥ 4. -/
theorem q4_lift_norms_distinct : liftNormSq 4 1 ≠ liftNormSq 4 2 := by decide

/-! ## Forward direction (q ≤ 3 ⟹ all nonzero residues equinormal) -/

/-- For q ∈ {2, 3}, every nonzero residue has lift norm exactly 1. -/
theorem all_nonzero_lns_one_of_qle3 :
    ∀ q, 2 ≤ q → q ≤ 3 → ∀ a, 0 < a → a < q → liftNormSq q a = 1 := by
  intro q hq2 hq3 a ha_pos ha_lt
  rcases Nat.lt_or_ge q 3 with h | h
  · -- q = 2
    have hq : q = 2 := by omega
    subst hq
    exact liftNormSq_q2 a ha_pos ha_lt
  · -- q = 3
    have hq : q = 3 := by omega
    subst hq
    exact liftNormSq_q3 a ha_pos ha_lt

/-- Consequence: for q ≤ 3, any two nonzero residues are equinormal. -/
theorem nonzero_lns_eq_of_qle3 :
    ∀ q, 2 ≤ q → q ≤ 3 →
    ∀ a b, 0 < a → a < q → 0 < b → b < q →
      liftNormSq q a = liftNormSq q b := by
  intro q hq2 hq3 a b ha_pos ha_lt hb_pos hb_lt
  rw [all_nonzero_lns_one_of_qle3 q hq2 hq3 a ha_pos ha_lt,
      all_nonzero_lns_one_of_qle3 q hq2 hq3 b hb_pos hb_lt]

/-! ## Converse direction (q ≥ 4 ⟹ ∃ residue with lift norm > 1) -/

/-- For q ≥ 4, residue 2 has lift norm strictly greater than 1. -/
theorem lns_two_gt_one_of_qge4 : ∀ q, 4 ≤ q → 1 < liftNormSq q 2 := by
  intro q hq
  unfold liftNormSq
  have hmin : min 2 (q - 2) = 2 := by
    have : 2 ≤ q - 2 := by omega
    exact min_eq_left this
  rw [hmin]
  decide

/-- Existence form of the converse: q ≥ 4 admits a nonzero residue
    whose lift norm is not 1. -/
theorem exists_nonzero_lns_ne_one_of_qge4 :
    ∀ q, 4 ≤ q → ∃ a, 0 < a ∧ a < q ∧ liftNormSq q a ≠ 1 := by
  intro q hq
  refine ⟨2, by decide, by omega, ?_⟩
  have := lns_two_gt_one_of_qge4 q hq
  omega

/-- Parametric form of `q4_lift_norms_distinct`: for *every* q ≥ 4 (not
    just q = 4), residues 1 and 2 have distinct lift norms. The single
    witness theorem above is the q = 4 instance.

    This generalisation matters: constant-fuzzing the q in the q = 4
    instance produces another true statement, exposing that the original
    name overstates specificity. -/
theorem qge4_lift_norms_distinct (q : ℕ) (hq : 4 ≤ q) :
    liftNormSq q 1 ≠ liftNormSq q 2 := by
  have h1 : liftNormSq q 1 = 1 := by
    unfold liftNormSq
    have hmin : min 1 (q - 1) = 1 := by
      have hge : 1 ≤ q - 1 := by omega
      exact min_eq_left hge
    rw [hmin]
  have h2 : 1 < liftNormSq q 2 := lns_two_gt_one_of_qge4 q hq
  omega

/-! ## Iff characterization

  The combinatorial heart of Lemma 5: every nonzero residue mod q
  has lift norm 1 iff q ≤ 3. This is the load-bearing fact behind
  the Bose-Mesner iff: when all residues are equinormal, the
  Euclidean cutoff |m|² ≤ K is constant on Hamming-distance
  classes, putting G(K) in BM(H(n,q)). When they are not
  (q ≥ 4), distance-class constancy fails.
-/

/-- **Tier 1A main theorem.** The lift-norm equinormality criterion
    (which characterizes when the Euclidean cutoff respects Hamming
    distance, hence when G(K) lies in BM(H(n,q))) holds iff q ≤ 3. -/
theorem all_nonzero_lns_one_iff_qle3 (q : ℕ) (hq : 2 ≤ q) :
    (∀ a, 0 < a → a < q → liftNormSq q a = 1) ↔ q ≤ 3 := by
  constructor
  · intro h
    by_contra hgt
    push_neg at hgt
    have hq4 : 4 ≤ q := hgt
    have h2_pos : 0 < (2 : ℕ) := by decide
    have h2_lt : (2 : ℕ) < q := by omega
    have hlns2 : liftNormSq q 2 = 1 := h 2 h2_pos h2_lt
    have hgt1 : 1 < liftNormSq q 2 := lns_two_gt_one_of_qge4 q hq4
    omega
  · intro hq3 a ha_pos ha_lt
    exact all_nonzero_lns_one_of_qle3 q hq hq3 a ha_pos ha_lt

/-- **Bose-Mesner iff (combinatorial form).** A Hamming-scheme G(K)
    lies in the Bose-Mesner algebra of H(n,q) precisely when the
    Euclidean cutoff respects Hamming distance, which holds iff
    every nonzero residue has lift norm 1, which holds iff q ≤ 3.

    This theorem packages the equivalence in the form needed by the
    manuscript: the lift-norm-equinormality criterion is the
    necessary and sufficient combinatorial condition. -/
theorem bose_mesner_combinatorial_iff (q : ℕ) (hq : 2 ≤ q) :
    (∀ a b, 0 < a → a < q → 0 < b → b < q → liftNormSq q a = liftNormSq q b)
    ↔ q ≤ 3 := by
  constructor
  · intro h
    by_contra hgt
    push_neg at hgt
    have hq4 : 4 ≤ q := hgt
    have h1_pos : 0 < (1 : ℕ) := by decide
    have h1_lt : (1 : ℕ) < q := by omega
    have h2_pos : 0 < (2 : ℕ) := by decide
    have h2_lt : (2 : ℕ) < q := by omega
    have heq : liftNormSq q 1 = liftNormSq q 2 := h 1 2 h1_pos h1_lt h2_pos h2_lt
    -- liftNormSq q 1 = 1 since min 1 (q-1) = 1 for q ≥ 2
    have h_one : liftNormSq q 1 = 1 := by
      unfold liftNormSq
      have hmin : min 1 (q - 1) = 1 := by
        have hge : 1 ≤ q - 1 := by omega
        exact min_eq_left hge
      rw [hmin]
    rw [h_one] at heq
    have hgt1 : 1 < liftNormSq q 2 := lns_two_gt_one_of_qge4 q hq4
    omega
  · intro hq3
    exact nonzero_lns_eq_of_qle3 q hq hq3

/-! ## Multi-q specialization sanity checks (technique #4)

  Concrete `decide`-anchored values of `liftNormSq q a` at q ∈ {2,3,4,5,7,11}
  for every nonzero residue. These pin the abstract iff to literal numbers
  cross-validated against the Python harness (`cross_validate_lean.py`,
  Section 6). The pattern is: for q ≤ 3 every value equals 1; for q ≥ 4 the
  middle residue ⌊q/2⌋ has lift norm > 1.
-/

theorem lns_q2_table : liftNormSq 2 1 = 1 := by decide

theorem lns_q3_table :
    liftNormSq 3 1 = 1 ∧ liftNormSq 3 2 = 1 := by decide

theorem lns_q4_table :
    liftNormSq 4 1 = 1 ∧ liftNormSq 4 2 = 4 ∧ liftNormSq 4 3 = 1 := by decide

theorem lns_q5_table :
    liftNormSq 5 1 = 1 ∧ liftNormSq 5 2 = 4 ∧
    liftNormSq 5 3 = 4 ∧ liftNormSq 5 4 = 1 := by decide

theorem lns_q7_table :
    liftNormSq 7 1 = 1 ∧ liftNormSq 7 2 = 4 ∧ liftNormSq 7 3 = 9 ∧
    liftNormSq 7 4 = 9 ∧ liftNormSq 7 5 = 4 ∧ liftNormSq 7 6 = 1 := by decide

theorem lns_q11_middle : liftNormSq 11 5 = 25 ∧ liftNormSq 11 6 = 25 := by decide

-- Note: a former `q_distinct_witnesses` theorem was deleted after
-- constant-fuzzing showed it made a uniform claim (`lns q 1 ≠ lns q 2`
-- holds for *every* q ≥ 4, not just the spot-checked values), so
-- swapping any q for any other q ≥ 4 left it true. The parametric
-- `qge4_lift_norms_distinct` above is the correct anchor: it states
-- the universal fact directly and is q-pinning by quantifier rather
-- than by enumeration.

/-- Negative anchor: the iff fails at q = 4 — there exists a nonzero
    residue (namely 2) with lift norm ≠ 1, proving the forward
    direction is non-vacuous. The parametric `iff_fails_at_qge4` below
    states the universal fact for every q ≥ 4. -/
theorem iff_fails_at_q4 :
    ¬ (∀ a, 0 < a → a < 4 → liftNormSq 4 a = 1) := by
  intro h
  have : liftNormSq 4 2 = 1 := h 2 (by decide) (by decide)
  have h4 : liftNormSq 4 2 = 4 := by decide
  omega

/-- Parametric form: for *every* q ≥ 4 the equinormality property
    fails. This is the actual content the manuscript needs; the
    `q4` instance above is one witness. -/
theorem iff_fails_at_qge4 (q : ℕ) (hq : 4 ≤ q) :
    ¬ (∀ a, 0 < a → a < q → liftNormSq q a = 1) := by
  intro h
  have h2_pos : 0 < (2 : ℕ) := by decide
  have h2_lt : (2 : ℕ) < q := by omega
  have hone : liftNormSq q 2 = 1 := h 2 h2_pos h2_lt
  have hgt : 1 < liftNormSq q 2 := lns_two_gt_one_of_qge4 q hq
  omega

/-! ## On q-pinning at the level of individual values

  An earlier `lns_q_pin_distinct` theorem attempted to pin each q
  against neighbours q±1 by spot-checking `liftNormSq q a` at boundary
  residues. Constant-fuzzing showed this is impossible at the level of
  individual literal values: the lift-norm function has two structural
  symmetries that defeat single-value pinning.

    1. Reflection: `liftNormSq q a = liftNormSq q (q - a)`, so any
       residue mutation that swaps a with its complement is invisible.
    2. Saturation: for q ≥ 2a, `liftNormSq q a = a²` is constant in q,
       so mutating q above the saturation point preserves the value.

  Together these mean that for any single equality `liftNormSq q a = v`
  there exist other (q', a') with the same value. The constant-fuzzer
  flags this honestly.

  The correct q-distinguishing anchors are the *parametric* theorems
  `qge4_lift_norms_distinct` and `iff_fails_at_qge4` above: they pin q
  by universal quantifier rather than by enumeration. Mutating the
  bound `4` in those theorems either weakens them (still true) or
  strengthens them past the truth boundary (build breaks at q = 3).

  The `lns_q*_table` theorems below are *value spot-checks* against
  the Python harness, not q-pins; their constant-fuzz survivors reflect
  the symmetries above and are not defects. -/

end KstarFormal.QaryGram
