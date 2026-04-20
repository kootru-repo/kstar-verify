/-
  KstarFormal.Combinatorics.GhzNonCoverage — GHZ stabilizer non-coverage by K*(4)
  K* Verification: Krawtchouk spectral correspondence

  Locks the load-bearing negative result behind the manuscript's GHZ
  reanalysis (Sec ghz, footnote $^\S$ at L1924): the K* operator set
  for n=4 omits exactly five of the GHZ stabilizers.

  The 5 missing operators are:
    XXYY, XYYX, YXXY, YYXX, ZZZZ

  The full GHZ stabilizer group has 16 nonzero Pauli operators
  (the projector |GHZ><GHZ| has these as nonzero coefficients).
  K*(4) covers 11 of them; the 5 above are omitted.

  Encoding: each Pauli is a 4-element list of {0,1,2,3} = {I,X,Y,Z}.
  The 137 K* labels are embedded verbatim from the certified output of
  `core.select_kstar_paulis(4)` in the upstream Python implementation.

  The proof is `decide` against the embedded list. This converts the
  manuscript's GHZ resolution from a script-verified claim into a
  Lean theorem that depends only on the embedded constants.

  Status: Tier 1 formalization, sorry-free.
-/
import Mathlib.Data.List.Basic
import Mathlib.Tactic.NormNum

namespace KstarFormal.GhzNonCoverage

/-! ## The 137 K* labels at n = 4 (encoded as [I=0, X=1, Y=2, Z=3]) -/

/-- The 137 K* Pauli labels for n=4, q=2 — verbatim from
    `core.select_kstar_paulis(4)`. Each entry is a 4-element list with
    values in {0,1,2,3} = {I,X,Y,Z}. -/
def kstar_labels_n4 : List (List Nat) := [
  [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 2], [0, 0, 0, 3],
  [0, 0, 1, 0], [0, 0, 1, 1], [0, 0, 1, 2], [0, 0, 1, 3],
  [0, 0, 2, 0], [0, 0, 2, 1], [0, 0, 2, 2], [0, 0, 2, 3],
  [0, 0, 3, 0], [0, 0, 3, 1], [0, 0, 3, 2], [0, 0, 3, 3],
  [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 0, 2], [0, 1, 0, 3],
  [0, 1, 1, 0], [0, 1, 1, 1], [0, 1, 1, 2], [0, 1, 1, 3],
  [0, 1, 2, 0], [0, 1, 2, 1], [0, 1, 2, 2], [0, 1, 2, 3],
  [0, 1, 3, 0], [0, 1, 3, 1], [0, 1, 3, 2], [0, 1, 3, 3],
  [0, 2, 0, 0], [0, 2, 0, 1], [0, 2, 0, 2], [0, 2, 0, 3],
  [0, 2, 1, 0], [0, 2, 1, 1], [0, 2, 1, 2], [0, 2, 1, 3],
  [0, 2, 2, 0], [0, 2, 2, 1], [0, 2, 2, 2], [0, 2, 2, 3],
  [0, 2, 3, 0], [0, 2, 3, 1], [0, 2, 3, 2], [0, 2, 3, 3],
  [0, 3, 0, 0], [0, 3, 0, 1], [0, 3, 0, 2], [0, 3, 0, 3],
  [0, 3, 1, 0], [0, 3, 1, 1], [0, 3, 1, 2], [0, 3, 1, 3],
  [0, 3, 2, 0], [0, 3, 2, 1], [0, 3, 2, 2], [0, 3, 2, 3],
  [0, 3, 3, 0], [0, 3, 3, 1], [0, 3, 3, 2], [0, 3, 3, 3],
  [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 0, 2], [1, 0, 0, 3],
  [1, 0, 1, 0], [1, 0, 1, 1], [1, 0, 1, 2], [1, 0, 1, 3],
  [1, 0, 2, 0], [1, 0, 2, 1], [1, 0, 2, 2], [1, 0, 2, 3],
  [1, 0, 3, 0], [1, 0, 3, 1], [1, 1, 0, 0], [1, 1, 0, 1],
  [1, 1, 1, 0], [1, 1, 1, 1], [1, 1, 2, 3], [1, 2, 0, 0],
  [1, 2, 0, 1], [1, 2, 1, 0], [1, 2, 1, 2], [1, 2, 3, 1],
  [1, 3, 0, 0], [1, 3, 0, 1], [1, 3, 1, 0], [1, 3, 1, 3],
  [1, 3, 3, 2], [2, 0, 0, 0], [2, 0, 0, 1], [2, 0, 0, 2],
  [2, 0, 0, 3], [2, 0, 1, 0], [2, 0, 1, 1], [2, 0, 2, 0],
  [2, 0, 2, 1], [2, 0, 3, 0], [2, 0, 3, 1], [2, 1, 0, 0],
  [2, 1, 0, 1], [2, 1, 1, 0], [2, 1, 2, 1], [2, 1, 3, 3],
  [2, 2, 0, 0], [2, 2, 0, 1], [2, 2, 1, 0], [2, 2, 2, 2],
  [2, 3, 0, 0], [2, 3, 0, 1], [2, 3, 1, 0], [2, 3, 1, 1],
  [2, 3, 2, 3], [3, 0, 0, 0], [3, 0, 0, 1], [3, 0, 0, 2],
  [3, 0, 0, 3], [3, 0, 1, 0], [3, 0, 1, 1], [3, 0, 2, 0],
  [3, 0, 2, 1], [3, 0, 3, 0], [3, 0, 3, 1], [3, 1, 0, 0],
  [3, 1, 0, 1], [3, 1, 1, 0], [3, 1, 1, 2], [3, 1, 3, 1],
  [3, 2, 0, 0], [3, 2, 1, 3], [3, 2, 3, 2], [3, 3, 0, 0],
  [3, 3, 2, 1]
]

/-! ## The 5 GHZ stabilizers absent from K*(4) -/

/-- XXYY = [1, 1, 2, 2] -/
def ghz_XXYY : List Nat := [1, 1, 2, 2]
/-- XYYX = [1, 2, 2, 1] -/
def ghz_XYYX : List Nat := [1, 2, 2, 1]
/-- YXXY = [2, 1, 1, 2] -/
def ghz_YXXY : List Nat := [2, 1, 1, 2]
/-- YYXX = [2, 2, 1, 1] -/
def ghz_YYXX : List Nat := [2, 2, 1, 1]
/-- ZZZZ = [3, 3, 3, 3] -/
def ghz_ZZZZ : List Nat := [3, 3, 3, 3]

/-- The 5 GHZ stabilizers absent from K*(4). -/
def ghz_missing : List (List Nat) :=
  [ghz_XXYY, ghz_XYYX, ghz_YXXY, ghz_YYXX, ghz_ZZZZ]

/-! ## Sanity anchors on the embedded list -/

/-- The K* list has exactly 137 entries. -/
theorem kstar_count_n4 : kstar_labels_n4.length = 137 := by native_decide

/-- Every K* label has length 4 (the n=4 condition). -/
theorem kstar_all_length_four :
    kstar_labels_n4.all (fun l => l.length = 4) = true := by native_decide

/-- Every entry of every K* label is in {0,1,2,3}. -/
theorem kstar_all_alphabet :
    kstar_labels_n4.all (fun l => l.all (fun x => x < 4)) = true := by native_decide

/-- The 5 missing-stabilizer list has exactly 5 entries. -/
theorem ghz_missing_count : ghz_missing.length = 5 := by native_decide

/-! ## The non-coverage theorem -/

/-- **GHZ stabilizer non-coverage (load-bearing for the GHZ resolution).**
    None of the five GHZ stabilizers XXYY, XYYX, YXXY, YYXX, ZZZZ
    appear in the K*(4) operator set. -/
theorem ghz_missing_disjoint_from_kstar :
    ghz_missing.all (fun l => ¬ l ∈ kstar_labels_n4) = true := by native_decide

/-- Equivalent statement form: the intersection is empty. -/
theorem ghz_missing_kstar_intersection_empty :
    (ghz_missing.filter (· ∈ kstar_labels_n4)) = [] := by native_decide

/-- Per-operator non-membership, individually checkable. -/
theorem ghz_XXYY_not_in_kstar : ¬ ghz_XXYY ∈ kstar_labels_n4 := by native_decide
theorem ghz_XYYX_not_in_kstar : ¬ ghz_XYYX ∈ kstar_labels_n4 := by native_decide
theorem ghz_YXXY_not_in_kstar : ¬ ghz_YXXY ∈ kstar_labels_n4 := by native_decide
theorem ghz_YYXX_not_in_kstar : ¬ ghz_YYXX ∈ kstar_labels_n4 := by native_decide
theorem ghz_ZZZZ_not_in_kstar : ¬ ghz_ZZZZ ∈ kstar_labels_n4 := by native_decide

end KstarFormal.GhzNonCoverage
