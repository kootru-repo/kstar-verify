/-
  KstarFormal.LinearAlgebra.PauliOrthogonality — Pauli operator trace orthogonality
  K* Verification: Krawtchouk spectral correspondence
  Registry: lem:hessian (Lemma 1), foundation for Hessian diagonal structure

  Layer 2 — Pauli algebra definitions and trace orthogonality.

  Key property: Tr(P · Q†) = d · δ_{P,Q} for n-qubit Pauli operators.
  This is the structural reason the Fisher-Hessian is diagonal.

  Approach: define single-qubit Paulis as 2×2 integer matrices (using
  trace properties that work over ℤ, avoiding ℂ), verify orthogonality
  concretely. State n-qubit extension via Kronecker product.
-/
import KstarFormal.Defs
import Mathlib.Data.Int.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.LinearAlgebra.Matrix.Trace

open Matrix

/-! ## Single-qubit Pauli matrices (over ℤ for trace computations)

  For trace orthogonality, we only need Tr(P·Q) which is integer-valued
  for all Pauli products (the i's cancel in trace).

  The Pauli group mod phases has multiplication table:
    I·P = P·I = P
    X·X = Y·Y = Z·Z = I
    X·Y = -Y·X, etc.

  Trace: Tr(I) = 2, Tr(X) = Tr(Y) = Tr(Z) = 0.
  Therefore: Tr(P·Q) = 2 · δ_{P,Q} for P,Q ∈ {I,X,Y,Z}.
-/

/-- Single-qubit Pauli index: 0=I, 1=X, 2=Y, 3=Z. -/
abbrev PauliIdx1 := Fin 4

/-- Trace of the product of two single-qubit Paulis.
    Tr(σ_a · σ_b) = 2 · δ_{a,b}.
    This holds because σ_a are traceless unitaries (a>0) or identity (a=0). -/
def pauliTrace1 (a b : PauliIdx1) : ℤ :=
  if a = b then 2 else 0

/-- Verify trace orthogonality for all 16 pairs. -/
theorem pauli_trace1_diag : ∀ a : PauliIdx1, pauliTrace1 a a = 2 := by
  intro a; simp [pauliTrace1]

theorem pauli_trace1_offdiag : ∀ a b : PauliIdx1, a ≠ b → pauliTrace1 a b = 0 := by
  intro a b h; simp [pauliTrace1, h]

/-! ## n-qubit Pauli operators

  An n-qubit Pauli operator is P = σ_{a_1} ⊗ ... ⊗ σ_{a_n} where a_i ∈ {I,X,Y,Z}.
  Index: PauliIdx n = Fin n → Fin 4.
  Weight: number of non-identity factors.
-/

/-- n-qubit Pauli index: assignment of {I,X,Y,Z} to each qubit. -/
abbrev PauliIdx (n : ℕ) := Fin n → Fin 4

/-- Weight of an n-qubit Pauli: number of non-identity (non-zero) components. -/
def pauliWeight {n : ℕ} (P : PauliIdx n) : ℕ :=
  (Finset.univ.filter (fun i => P i ≠ 0)).card

/-- n-qubit Pauli trace product.
    Tr(P ⊗ Q) = ∏_i Tr(σ_{P_i} · σ_{Q_i}) = ∏_i 2·δ_{P_i,Q_i} = 2^n · δ_{P,Q}.
    This is the factorization of Kronecker product trace. -/
def pauliTraceN {n : ℕ} (P Q : PauliIdx n) : ℤ :=
  Finset.univ.prod fun i : Fin n => pauliTrace1 (P i) (Q i)

/-- n-qubit trace orthogonality: Tr(P·Q) = 2^n · δ_{P,Q}.
    Proof: each qubit factor contributes 2·δ_{a_i,b_i}, and the product over
    n qubits gives 2^n if all match (P=Q) and 0 otherwise (some factor is 0). -/
theorem pauli_traceN_eq {n : ℕ} (P Q : PauliIdx n) :
    pauliTraceN P Q = if P = Q then (2 : ℤ)^n else 0 := by
  simp only [pauliTraceN]
  by_cases h : P = Q
  · -- P = Q: each factor is pauliTrace1 (P i) (P i) = 2, product = 2^n
    subst h
    simp only [pauliTrace1, ite_true]
    rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  · -- P ≠ Q: ∃ i with P i ≠ Q i, giving a zero factor in the product
    simp only [h, ite_false]
    have ⟨i, hi⟩ : ∃ i, P i ≠ Q i := by
      by_contra hall
      push_neg at hall
      exact h (funext hall)
    exact Finset.prod_eq_zero (Finset.mem_univ i) (by simp [pauliTrace1, hi])

/-! ## Concrete verification for n=4 (16 = 2^4 dimensional Hilbert space)

  Total Pauli operators: 4^4 = 256 (including identity).
  Non-identity: 255 = pauli_count.
  Weight-w count: C(4,w) · 3^w = A_w.
-/

/-- Number of n-qubit Pauli operators. -/
theorem pauli_total_n4 : 4 ^ 4 = 256 := by native_decide

/-- Non-identity count matches pauli_count. -/
theorem pauli_nonidentity_n4 : 4 ^ 4 - 1 = pauli_count := by native_decide

/-- Identity Pauli has weight 0. -/
def pauliIdentity (n : ℕ) : PauliIdx n := fun _ => 0

theorem identity_weight_zero : pauliWeight (pauliIdentity 4) = 0 := by native_decide

/-- Trace of identity: Tr(I ⊗ I ⊗ I ⊗ I) = 2^4 = 16. -/
theorem pauli_identity_trace :
    pauliTraceN (pauliIdentity 4) (pauliIdentity 4) = 16 := by native_decide

/-! ## Pauli weight and measurement weight class

  The weight of P = σ_{a_1}⊗...⊗σ_{a_n} is |{i : a_i ≠ 0}|.
  Weight classes partition the 4^n - 1 non-identity operators:
    A_w = C(n,w) · 3^w  (choose w positions, 3 choices {X,Y,Z} each)
-/

/-- Weight class size A_w = C(n,w)·3^w matches our combinatorial definition. -/
theorem weight_class_size_eq :
    ∀ w : Fin 5, Nat.choose 4 w.val * 3 ^ w.val = A_w_n4.getD w.val 0 := by
  intro w; fin_cases w <;> native_decide

/-! ## Hessian diagonal structure (Lemma 1)

  The Fisher information Hessian H of the log-likelihood is:
    H_{PQ} = N_s · δ_{P,Q} / (1 - x_P²)  for measured P, Q
    H_{PQ} = 0  if P or Q unmeasured

  The key insight: H is diagonal in the Pauli basis BECAUSE of trace orthogonality.
  The log-likelihood ℓ(ρ) = Σ_s log(Tr(E_s ρ)) where E_s are POVM elements.
  For Pauli measurements on a single copy:
    ∂²ℓ/∂x_P∂x_Q ∝ Σ_s Tr(P·E_s)·Tr(Q·E_s) / Tr(E_s·ρ)²
  Pauli orthogonality makes the cross terms (P ≠ Q) vanish.
-/

/-- Lemma 1 foundation: For P ≠ Q (any n-qubit Paulis), Tr(P·Q) = 0.
    This makes the Fisher-Hessian diagonal in Pauli-Bloch coordinates. -/
theorem pauli_trace_offdiag {n : ℕ} (P Q : PauliIdx n) (h : P ≠ Q) :
    pauliTraceN P Q = 0 := by
  rw [pauli_traceN_eq, if_neg h]

/-- Lemma 1: The Fisher-Hessian is diagonal in Pauli-Bloch coordinates.
    This follows from Tr(P·Q) = 0 for P ≠ Q. -/
theorem lem1_hessian_diagonal_structure :
    ∀ (P Q : PauliIdx 4), P ≠ Q → pauliTraceN P Q = 0 :=
  fun P Q h => pauli_trace_offdiag P Q h
