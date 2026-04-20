import Lake
open Lake DSL

package KstarFormal where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

-- doc-gen4 must match the Lean toolchain (v4.29.0-rc8). Require it
-- BEFORE mathlib so mathlib's transitive pins take precedence for
-- shared dependencies (plausible, aesop, etc.) -- otherwise the
-- `lake exe cache get` hash check rejects mismatched revisions.
require «doc-gen4» from git
  "https://github.com/leanprover/doc-gen4" @ "v4.29.0-rc8"

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "master"

@[default_target]
lean_lib KstarFormal where
  srcDir := "."
