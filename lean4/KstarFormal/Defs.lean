/-
  KstarFormal.Defs — Shared definitions for K* formal verification
  K* Verification: Krawtchouk spectral correspondence
  Registry: shared definitions used across all layers
-/
import Mathlib.Data.Int.Basic
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Data.List.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum

open Finset BigOperators

/-! ## Constants for the n=4, q=2 Hamming scheme -/

/-- Number of qubits. -/
abbrev n_qubits : ℕ := 4

/-- Local dimension q (binary). -/
abbrev q_local : ℕ := 2

/-- Hilbert space dimension d = q^n = 16. -/
abbrev hilbert_dim : ℕ := 2 ^ n_qubits

/-- Total non-identity Pauli count: q^{2n} - 1 = 255. -/
abbrev pauli_count : ℕ := 2 ^ (2 * n_qubits) - 1

/-! ## Weight class sizes A_w = C(n,w) · (q²-1)^w -/

/-- Size of Pauli weight class w for n qubits, alphabet q. -/
def weightClassSize (n w : ℕ) (q : ℕ := 2) : ℕ :=
  Nat.choose n w * (q ^ 2 - 1) ^ w

/-- A_w for n=4, q=2: [1, 12, 54, 108, 81]. -/
def A_w_n4 : List ℕ := [1, 12, 54, 108, 81]

theorem A_w_n4_correct :
    A_w_n4 = [weightClassSize 4 0, weightClassSize 4 1,
              weightClassSize 4 2, weightClassSize 4 3,
              weightClassSize 4 4] := by
  native_decide

theorem A_w_n4_sum : A_w_n4.sum = 256 := by native_decide

/-! ## Parity weight of a Z^n vector -/

/-- Parity weight: number of odd entries. -/
def parityWeight {n : ℕ} (v : Fin n → ℤ) : ℕ :=
  (Finset.univ.filter (fun i => v i % 2 ≠ 0)).card
