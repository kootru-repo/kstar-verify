/-
  KstarFormal.Combinatorics.PatchDecomp — Phase G: patch decomposition
  K* Verification: Krawtchouk spectral correspondence

  Phase G (universal-n airtightness plan, 2026-04-17):
  formalize the compositional-architecture patch decomposition
  used by the paper's Sec. "8qubit" to extend K* = 5 saturation
  from monolithic n = 4 to larger systems n >= 4.

  Paper claim: "Cover the n qubits with a collection of 4-qubit
  patches; each patch saturates internally at K* = 5; the
  aggregate covers all weight-1 and weight-2 supports."

  Core content:
  * `Patch n` = 4-subset of `Fin n`                                (subtype)
  * `PatchCover n` = List (Patch n)                                (finite list)
  * `pauliSupport P` = Finset of indices where P is non-identity   (Finset (Fin n))
  * `coversPauli cover P` = P.support is contained in some patch   (Prop)
  * `patch_superset_of_small` : any |S| ≤ 4 extends to a 4-subset  (key lemma)

  Gate G (this file + UniversalKstarCompositional.lean):
  0 sorry, 0 new propositional axioms, universal over n >= 4.
-/
import KstarFormal.LinearAlgebra.PauliOrthogonality
import Mathlib.Data.Finset.Card

open Finset

/-! ## Pauli support as a Finset -/

/-- Support of an n-qubit Pauli: the set of indices with non-identity factor.
    By construction, `(pauliSupport P).card = pauliWeight P`. -/
def pauliSupport {n : ℕ} (P : PauliIdx n) : Finset (Fin n) :=
  Finset.univ.filter (fun i => P i ≠ 0)

theorem pauliSupport_card_eq_weight {n : ℕ} (P : PauliIdx n) :
    (pauliSupport P).card = pauliWeight P := rfl

/-! ## Patches: 4-qubit subsets of Fin n

  A patch is a 4-element subset of the index set Fin n.  The paper's
  compositional architecture uses a cover by patches so that every
  low-weight Pauli has its support contained in some patch, allowing
  per-patch K* = 5 saturation to compose.
-/

/-- A patch is a 4-element subset of `Fin n`. -/
def Patch (n : ℕ) : Type := {S : Finset (Fin n) // S.card = 4}

/-- A patch cover is a list of patches. -/
abbrev PatchCover (n : ℕ) : Type := List (Patch n)

/-- Predicate: the Pauli `P` is covered by the patch cover `cover`
    if its support is contained in at least one patch. -/
def coversPauli {n : ℕ} (cover : PatchCover n) (P : PauliIdx n) : Prop :=
  ∃ p ∈ cover, pauliSupport P ⊆ p.val

/-! ## Key extension lemma: any small support extends to a 4-patch

  When n >= 4, any `S : Finset (Fin n)` with `S.card <= 4` embeds
  into some 4-subset of `Fin n`.  This is the structural reason why
  a patch cover of Fin n can simultaneously cover every weight-≤-4
  Pauli support.
-/

/-- For n >= 4, every subset of Fin n with card ≤ 4 extends to a
    4-subset of Fin n.

    Proof: pick (4 - |S|) elements from the complement `univ \ S`
    using `Finset.exists_subset_card_eq`; their union with S has
    cardinality 4 via the disjoint-union card formula. -/
theorem patch_superset_of_small {n : ℕ} (h : 4 ≤ n)
    {S : Finset (Fin n)} (hS : S.card ≤ 4) :
    ∃ T : Finset (Fin n), S ⊆ T ∧ T.card = 4 := by
  have hcomp_card :
      ((Finset.univ : Finset (Fin n)) \ S).card = n - S.card := by
    rw [Finset.card_sdiff_of_subset (Finset.subset_univ S)]
    simp [Finset.card_univ, Fintype.card_fin]
  have hbound : 4 - S.card ≤ ((Finset.univ : Finset (Fin n)) \ S).card := by
    rw [hcomp_card]; omega
  obtain ⟨U, hUsub, hUcard⟩ := Finset.exists_subset_card_eq hbound
  refine ⟨S ∪ U, Finset.subset_union_left, ?_⟩
  have hdisj : Disjoint S U :=
    Finset.disjoint_of_subset_right hUsub Finset.disjoint_sdiff
  rw [Finset.card_union_of_disjoint hdisj, hUcard]
  omega

/-- Convenience wrapper: a Pauli with weight ≤ 2 has its support
    embeddable in a 4-patch whenever n >= 4. -/
theorem pauli_weight_le_two_extends_to_patch {n : ℕ} (h : 4 ≤ n)
    (P : PauliIdx n) (hw : pauliWeight P ≤ 2) :
    ∃ T : Finset (Fin n), pauliSupport P ⊆ T ∧ T.card = 4 := by
  have : (pauliSupport P).card ≤ 4 := by
    rw [pauliSupport_card_eq_weight]; omega
  exact patch_superset_of_small h this

/-- Convenience wrapper: a Pauli with weight ≤ 4 has its support
    embeddable in a 4-patch whenever n >= 4. -/
theorem pauli_weight_le_four_extends_to_patch {n : ℕ} (h : 4 ≤ n)
    (P : PauliIdx n) (hw : pauliWeight P ≤ 4) :
    ∃ T : Finset (Fin n), pauliSupport P ⊆ T ∧ T.card = 4 := by
  have : (pauliSupport P).card ≤ 4 := by
    rw [pauliSupport_card_eq_weight]; exact hw
  exact patch_superset_of_small h this
