# K* Formal Verification (Lean4)

Machine-verified proofs for the companion manuscript:
**"Krawtchouk spectral correspondence between Hamming schemes and quantum measurement design"**

## Status

| Layer | Scope | Status | Sorry count |
|-------|-------|--------|-------------|
| **1** | Pure combinatorics (K\*, eigenvalues, allocation, n=4 concrete) | Complete | 0 |
| **2** | Quantum linear algebra (Hessian, purity, fidelity, PSD contraction) | Complete (from axioms) | 0 |
| **3** | Probability (hypergeometric, Hoeffding) | Complete | 0 |

All three layers are sorry-free. The physical conclusions (F=0 witness pair, min F ≤ 1/2, HS bound) are proved from 4 clearly-labeled propositional axioms in `KstarFormal/Axioms.lean` — standard textbook results (Uhlmann, Schatten, fidelity-sum-complementary, fidelity-nonneg) that Mathlib does not yet expose. `#print axioms <theorem>` reveals the per-theorem axiom dependencies.

## Machine-Verified Claims (Layer 1)

- N₄(5) = 137 (lattice point count)
- c_w(K=5) = [9, 56, 24, 32, 16] (parity weight distribution)
- λ_w(K=5) = [144, 224, 64, 128, 256] (Gram eigenvalues)
- All λ_w > 0 (positive definiteness)
- κ(G) = 256/64 = 4 (condition number)
- M_w(K=5) = [1, 12, 54, 54, 16] (greedy redistribution)
- w_sat = 2 (saturation index)
- K\* = 5 (minimal cutoff)
- |S| ≥ 66 (Corollary 2 operator lower bound)
- λ_w(K=4) ≤ λ_w(K=5) (Lemma 3 monotonicity)
- Krawtchouk orthogonality for n=4

## Lean4 ↔ Paper Mapping

| Registry ID | Paper | Lean file | Status |
|-------------|-------|-----------|--------|
| `prop:spectral_q_main` | Lemma 5 | `LinearAlgebra/SpectralDecomp.lean` | Layer 1 |
| `app:completeness` | Lemma 6 | `LinearAlgebra/Eigenvalues.lean` | Layer 1 |
| `lem:monotone` | Lemma 3 | `LinearAlgebra/Monotonicity.lean` | Layer 1 |
| `cor:lower_bound` | Corollary 2 | `Combinatorics/WeightSaturation.lean` | Layer 1 |
| `lem:hessian` | Lemma 1 | `LinearAlgebra/PauliOrthogonality.lean` | Layer 2 (0 sorry) |
| `prop:purity_main` | Lemma 2 | `Quantum/PurityBound.lean` | Layer 2 (0 sorry) |
| `thm:basin` | Theorem 1 | `Quantum/BasinSeparation.lean` | Layer 2 (0 sorry) |
| `thm:spectral-char` | Theorem 2 | `Quantum/FidelityDichotomy.lean` | Layer 2 (0 sorry, from axioms) |
| `cor:approx_local` | Corollary 1 | `Quantum/BasinSeparation.lean` | Layer 2 (0 sorry) |
| `prop:coupon` | Lemma 4 | `Probability/Hypergeometric.lean` | Layer 3 (0 sorry) |
| `thm:asymptotic` | Theorem 3 | `Statements.lean` | Layer 2 (0 sorry) |

## Setup

### Prerequisites
1. Install [elan](https://github.com/leanprover/elan) (Lean version manager):
   ```powershell
   # PowerShell (admin)
   irm https://raw.githubusercontent.com/leanprover/elan/master/elan-init.ps1 | iex
   ```
2. Restart terminal to pick up PATH changes.
3. Install VS Code [Lean 4 extension](https://marketplace.visualstudio.com/items?itemName=leanprover.lean4).

### Build
```bash
cd lean4-formal
lake build        # First build fetches + compiles mathlib (~30-60 min)
```

### Audit
```bash
bash scripts/check_sorry.sh           # Report all sorry
bash scripts/check_sorry.sh --layer1  # Fail if Layer 1 has sorry
```

## Architecture

```
KstarFormal.lean              ← Root import (the certificate)
├── Defs.lean                  ← Shared definitions (n=4, q=2, A_w)
├── Combinatorics/
│   ├── Krawtchouk.lean        ← K_w(h;n) definition + n=4 matrix + orthogonality
│   ├── LatticeCount.lean      ← r₄(k), N₄(K)=137, c_w, shell decomposition
│   ├── GreedyRedist.lean      ← Greedy algorithm, M_w, monotonicity
│   └── WeightSaturation.lean  ← K*=5, w_sat=2, |S|≥66
├── LinearAlgebra/
│   ├── GramMatrix.lean        ← G(K), Hamming-constancy, distance values
│   ├── SpectralDecomp.lean    ← Eigenvalue formula, Krawtchouk transform
│   ├── Eigenvalues.lean       ← All-positive, rank=16, κ=4
│   ├── Monotonicity.lean      ← λ_w(4)≤λ_w(5), shell analysis
│   └── PauliOrthogonality.lean ← [Layer 2 stub] Tr(P·Q)=d·δ
├── Quantum/
│   ├── PurityBound.lean       ← [Layer 2 stub] Lemma 2
│   ├── FidelityDichotomy.lean ← [Layer 2 stub] Theorem 2(i),(iii)
│   └── BasinSeparation.lean   ← [Layer 2 stub] Theorem 1, Corollary 1
├── Probability/
│   └── Hypergeometric.lean    ← [Layer 3 stub] Lemma 4
├── Statements.lean            ← All 11 statements as Lean propositions
└── Certified.lean             ← Layer 1 assembly (sorry-free deliverable)
```

## Cross-References

- **proofs_registry.yaml**: `../external-proofs/proofs_registry.yaml` — dependency graph
- **certified_values.json**: `../sage-exact-verification/certified_values.json` — ground truth
- **Numerical verification**: `../verification/run_master.py` — 6-tier system (349 checks)
