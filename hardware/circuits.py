"""Circuit-building helpers for live-QPU K* replication.

`pauli_label_to_circuit(label, n_qubits, prep)` builds a Qiskit
circuit that (1) prepares a reference state, (2) rotates into the
requested Pauli measurement basis, (3) measures all qubits.

`ideal_state_for_prep(prep, n_qubits)` returns the noiseless
density matrix for the same reference state, used locally for
fidelity comparison against the reconstructed density matrix.

Keep this module dependency-light: qiskit imports live inside
`pauli_label_to_circuit` (the only function that builds circuits),
so the analyzer, the bounds loader, and tests that only exercise
`counts_to_expectation` or `ideal_state_for_prep` can import this
file without qiskit installed.
"""
from __future__ import annotations

import math

import numpy as np


def pauli_label_to_circuit(pauli_string: str, n_qubits: int, prep: str = "w"):
    """Build a QuantumCircuit that prepares `prep`, rotates into the
    measurement basis specified by `pauli_string`, and measures all
    qubits in the computational basis.

    pauli_string[i] in {'I','X','Y','Z'}; 'I' means no rotation (the
    identity contributes +1 to every shot, so we simply don't read
    that qubit's parity during post-processing).

    prep options:
      - 'w'     : 4-qubit W state via cascaded CRy/CNOT splitting
      - 'ghz'   : H on q0, CNOT chain q_i -> q_{i+1}
      - 'plus'  : H on every qubit (|+>^n)
      - 'bell'  : Bell pair on q0,q1; only valid for n_qubits >= 2
      - 'zero'  : |0>^n (no preparation gates)
    """
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

    # Name the classical register "meas" so SamplerV2 exposes the
    # counts as `pub_result.data.meas.get_counts()` -- the default
    # register name from `QuantumCircuit(n, n)` is "c", which would
    # force every caller to probe pub_result.data attributes by name.
    qr = QuantumRegister(n_qubits, "q")
    cr = ClassicalRegister(n_qubits, "meas")
    qc = QuantumCircuit(qr, cr)

    if prep == "w":
        if n_qubits != 4:
            raise ValueError("The W-state preparation shipped here is hard-"
                             "wired for n_qubits=4.  Extend the splitting "
                             "schedule for other sizes.")
        # Dicke D(4,1): start with |1000>, then cascade excitation across
        # q1/q2/q3 with splitting ratios 3/4, 2/3, 1/2.  Each split is a
        # controlled-Ry followed by a CNOT in the opposite direction.
        qc.x(0)
        qc.cry(2 * math.asin(math.sqrt(3 / 4)), 0, 1)
        qc.cx(1, 0)
        qc.cry(2 * math.asin(math.sqrt(2 / 3)), 1, 2)
        qc.cx(2, 1)
        qc.cry(2 * math.asin(math.sqrt(1 / 2)), 2, 3)
        qc.cx(3, 2)
    elif prep == "ghz":
        qc.h(0)
        for i in range(n_qubits - 1):
            qc.cx(i, i + 1)
    elif prep == "plus":
        for i in range(n_qubits):
            qc.h(i)
    elif prep == "bell":
        if n_qubits < 2:
            raise ValueError("bell requires n_qubits >= 2")
        qc.h(0)
        qc.cx(0, 1)
    elif prep == "zero":
        pass
    else:
        raise ValueError(f"unknown prep: {prep!r}")

    qc.barrier()

    for i, basis in enumerate(pauli_string):
        if basis == "X":
            qc.h(i)
        elif basis == "Y":
            qc.sdg(i)
            qc.h(i)
        elif basis in ("I", "Z"):
            pass
        else:
            raise ValueError(f"pauli_string contains non-Pauli char {basis!r}")

    qc.measure(range(n_qubits), range(n_qubits))
    return qc


def counts_to_expectation(counts: dict, pauli_string: str, n_qubits: int) -> float:
    """Convert raw counts to <P> for a single Pauli operator.

    Only qubits with pauli_string[i] != 'I' contribute to the parity;
    bitstrings are reversed once to undo Qiskit's little-endian
    convention so parity indexing matches pauli_string's qubit order.
    """
    total = sum(counts.values())
    if total == 0:
        return 0.0
    active_qubits = [i for i, p in enumerate(pauli_string) if p != "I"]
    if not active_qubits:
        return 1.0
    expectation = 0.0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        parity = sum(int(bits[q]) for q in active_qubits) % 2
        expectation += ((-1) ** parity) * count
    return expectation / total


def ideal_state_for_prep(prep: str, n_qubits: int) -> np.ndarray:
    """Return the noiseless density matrix for `prep`."""
    if prep == "w":
        psi = np.zeros(2 ** n_qubits, dtype=complex)
        for i in range(n_qubits):
            idx = 1 << (n_qubits - 1 - i)
            psi[idx] = 1.0 / math.sqrt(n_qubits)
    elif prep == "ghz":
        psi = np.zeros(2 ** n_qubits, dtype=complex)
        psi[0] = 1.0 / math.sqrt(2)
        psi[-1] = 1.0 / math.sqrt(2)
    elif prep == "plus":
        plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
        psi = plus
        for _ in range(n_qubits - 1):
            psi = np.kron(psi, plus)
    elif prep == "bell":
        psi = np.zeros(2 ** n_qubits, dtype=complex)
        psi[0] = 1.0 / math.sqrt(2)
        psi[3] = 1.0 / math.sqrt(2)
    elif prep == "zero":
        psi = np.zeros(2 ** n_qubits, dtype=complex)
        psi[0] = 1.0
    else:
        raise ValueError(f"unknown prep: {prep!r}")
    return np.outer(psi, psi.conj())
