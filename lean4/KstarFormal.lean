/-
  KstarFormal — Root import for K* formal verification
  K* Verification: Krawtchouk spectral correspondence

  This is the top-level entry point. Importing this file imports
  the entire formal verification project.

  If `lake build` succeeds on this file, all formal proofs compile.
  Use `scripts/check_sorry.sh` to audit for remaining sorry's.
-/

-- Layer 1: Pure combinatorics (sorry-free target)
import KstarFormal.Defs
import KstarFormal.Combinatorics.Krawtchouk
import KstarFormal.Combinatorics.LatticeCount
import KstarFormal.Combinatorics.LatticeCountGeneric
import KstarFormal.Combinatorics.CwClosedForm
import KstarFormal.Combinatorics.GreedyRedist
import KstarFormal.Combinatorics.PatchDecomp
import KstarFormal.Combinatorics.WeightSaturation
import KstarFormal.LinearAlgebra.GramMatrix
import KstarFormal.LinearAlgebra.SpectralDecomp
import KstarFormal.LinearAlgebra.Eigenvalues
import KstarFormal.LinearAlgebra.Monotonicity

import KstarFormal.Combinatorics.KstarLabels
import KstarFormal.Combinatorics.DelsarteDual
import KstarFormal.Combinatorics.KFullScaling

-- Layer 1 assembly
import KstarFormal.Certified

-- All 11 statements
import KstarFormal.Statements

-- Phase B (2026-04-17): universal-n / universal-K wrappers
import KstarFormal.Universal

-- Phase C (2026-04-17): universal thm:spectral-char parts (i)-(iii)
import KstarFormal.UniversalSpectralChar

-- Phase C part (iv) (2026-04-17): K*=5 saturation at n in {4..8}
import KstarFormal.UniversalKstar

-- Phase G (2026-04-17): compositional patch saturation for all n >= 4
import KstarFormal.UniversalKstarCompositional

-- Phase E (2026-04-17): adversarial refinement for Phases B/C
import KstarFormal.AdversarialPhaseE

-- Layer 2: Quantum linear algebra (proved from axioms)
import KstarFormal.LinearAlgebra.PauliOrthogonality
import KstarFormal.Quantum.PurityBound
import KstarFormal.Quantum.FidelityDichotomy
import KstarFormal.Quantum.PSDContraction
import KstarFormal.Quantum.BasinSeparation

-- Layer 3: Probability (proved)
import KstarFormal.Probability.Hypergeometric
import KstarFormal.Probability.Hoeffding

-- Adversarial stress tests (101 boundary/negative/cross-validation/hedged-envelope checks)
import KstarFormal.Adversarial

-- Computational recovery of all 5 axioms + 4 load-bearing wiring theorems
-- (Phase 4a, 2026-04-07: 3 inconsistent axioms removed.
--  v2 refinement, 2026-04-08: parseval_pauli_bloch promoted axiom→theorem.
--  See Axioms.lean header for the full inventory.)
import KstarFormal.AxiomFoundation
