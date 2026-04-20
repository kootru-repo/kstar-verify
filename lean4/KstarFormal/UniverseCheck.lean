/-
UniverseCheck.lean — Aspect 6 (universe polymorphism stress)
============================================================

Walks the entire Lean environment after importing all of KstarFormal,
extracts every constant whose name begins with "KstarFormal", and reports
its `levelParams` (universe parameters). PASS if every KstarFormal
constant has empty `levelParams` (i.e. universe-monomorphic).

Rationale: a universe-polymorphic theorem can have a hidden semantic
mismatch with its informal claim if the universe is silently specialized
at the use site. Our manuscript claims live entirely in `Type 0`
(`ℕ`, `ℤ`, `ℚ`, `Fin n`, `PauliIdx n`, plus the opaque
`DensityMatrix : Type`). This file makes that claim machine-checkable.

Run with:
  lake env lean KstarFormal/UniverseCheck.lean
-/
import KstarFormal

open Lean

/-- Auto-generated suffixes that Lean produces for inductives, structures,
    and pattern-matching defs. These are NOT user-written theorems and may
    inherit universe polymorphism mechanically (e.g. `match_1` returning into
    `Sort u_1`); excluding them isolates user-written declarations. -/
def autoGenSuffixes : List String :=
  [".rec", ".recOn", ".casesOn", ".brecOn", ".below", ".binductionOn",
   ".noConfusion", ".noConfusionType", ".ndrec", ".ndrecOn",
   ".injEq", ".mk.injEq", ".sizeOf_spec", ".mk.sizeOf_spec",
   ".eq_def", ".mk", ".mk.sizeOf", ".sizeOf"]

def isAutoGen (n : Name) : Bool :=
  let s := n.toString
  -- match_N equation lemmas
  if (s.splitOn ".match_").length ≥ 2 then true
  -- _eq_N / _proof_N
  else if (s.splitOn "._eq_").length ≥ 2 then true
  else if (s.splitOn "._proof_").length ≥ 2 then true
  -- standard inductive/structure auto-suffixes
  else autoGenSuffixes.any (fun suf =>
    let parts := s.splitOn suf
    parts.length ≥ 2 ∧ (parts[parts.length - 1]?).getD "_" == "")

/-- Modules to exclude from the check (this file itself, since our local
    `CheckState` structure produces auto-projections that would otherwise
    pollute the report). -/
def excludedModules : List String := ["KstarFormal.UniverseCheck"]

structure CheckState where
  total : Nat := 0
  sample : Array Name := #[]
  poly : Array (Name × List Name) := #[]

def runUniverseCheck : MetaM Unit := do
  let env ← getEnv
  let st := env.constants.fold (init := ({} : CheckState)) fun st name info =>
    if name.isInternal then st
    else if isAutoGen name then st
    else
      -- Filter by module of origin: keep only constants defined in a
      -- module whose name starts with "KstarFormal".
      match env.getModuleFor? name with
      | none => st
      | some modName =>
        let modStr := modName.toString
        if !((modStr.splitOn "KstarFormal").length ≥ 2) then st
        else if excludedModules.contains modStr then st
        else
          let st := { st with total := st.total + 1 }
          let st := if st.sample.size < 8 then { st with sample := st.sample.push name } else st
          if info.levelParams.isEmpty then st
          else { st with poly := st.poly.push (name, info.levelParams) }
  IO.println s!"[universe-check] sample names: {st.sample.toList}"
  let totalChecked := st.total
  let polymorphic := st.poly
  IO.println s!"[universe-check] scanned {totalChecked} KstarFormal constants"
  if polymorphic.isEmpty then
    IO.println "[universe-check] PASS: every KstarFormal constant is universe-monomorphic"
  else
    IO.println s!"[universe-check] FAIL: {polymorphic.size} polymorphic constants:"
    for (n, ls) in polymorphic do
      IO.println s!"  {n}  ::  {ls}"

#eval runUniverseCheck
