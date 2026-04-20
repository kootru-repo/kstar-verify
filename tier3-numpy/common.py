"""Shared numerical primitives (mpmath/numpy) for K* Verification verification.

Tier 2 uses mpmath for high-precision arithmetic and numpy for
matrix operations on density matrices and Pauli operators.
"""

import numpy as np
from itertools import product as cartesian


# ---- Pauli matrices ----

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]


def pauli_tensor(indices):
    """Build n-qubit Pauli operator from index tuple.

    indices: tuple of 0,1,2,3 -> I,X,Y,Z at each position.
    Returns 2^n x 2^n numpy array.
    """
    result = PAULIS[indices[0]]
    for idx in indices[1:]:
        result = np.kron(result, PAULIS[idx])
    return result


def pauli_weight(indices):
    """Hamming weight of a Pauli index tuple."""
    return sum(1 for x in indices if x != 0)


def all_paulis(n):
    """Generate all n-qubit Pauli index tuples."""
    return list(cartesian(range(4), repeat=n))


# ---- Standard quantum states ----

def w_state(n):
    """n-qubit W state: (|10...0> + |01...0> + ... + |0...01>) / sqrt(n).
    Returns density matrix.
    """
    dim = 2**n
    psi = np.zeros(dim, dtype=complex)
    for i in range(n):
        idx = 1 << (n - 1 - i)
        psi[idx] = 1.0 / np.sqrt(n)
    return np.outer(psi, psi.conj())


def ghz_state(n):
    """n-qubit GHZ state: (|0...0> + |1...1>) / sqrt(2).
    Returns density matrix.
    """
    dim = 2**n
    psi = np.zeros(dim, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    psi[dim - 1] = 1.0 / np.sqrt(2)
    return np.outer(psi, psi.conj())


def product_state_plus(n):
    """|+>^n state. Returns density matrix."""
    psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = psi_plus
    for _ in range(n - 1):
        psi = np.kron(psi, psi_plus)
    return np.outer(psi, psi.conj())


def one_local_state(n):
    """State with W(rho)=1: |0><0| on qubit 1, maximally mixed on rest.

    Only weight-1 Z_1 expectation is nonzero.
    """
    rho_0 = np.array([[1, 0], [0, 0]], dtype=complex)  # |0><0|
    rho_mixed = np.eye(2, dtype=complex) / 2
    rho = rho_0
    for _ in range(n - 1):
        rho = np.kron(rho, rho_mixed)
    return rho


def two_local_state(n):
    """State with W(rho)=2: Bell pair on qubits 1,2, maximally mixed on rest.

    Nonzero expectations: XX, YY, ZZ on qubits 1,2 (weight-2).
    """
    bell = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    rho_bell = np.outer(bell, bell.conj())
    rho_mixed = np.eye(2, dtype=complex) / 2
    rho = rho_bell
    for _ in range(n - 2):
        rho = np.kron(rho, rho_mixed)
    return rho


def depolarise(rho, p):
    """Apply depolarising noise: rho -> (1-p)*rho + p*I/d."""
    d = rho.shape[0]
    return (1 - p) * rho + p * np.eye(d) / d


# ---- Pauli expectation values ----

def pauli_expectations(rho, n):
    """Compute all Pauli expectation values tr(P*rho).

    Returns dict: pauli_tuple -> expectation value.
    """
    expectations = {}
    for indices in all_paulis(n):
        P = pauli_tensor(indices)
        val = np.real(np.trace(P @ rho))
        expectations[indices] = val
    return expectations


# ---- K* operator set ----

def kstar_operator_indices(n, K=5):
    """Return set of Pauli index tuples in the K*-structured set.

    Uses weight-class allocation: saturate weight 0, 1, 2 (all ops),
    then fill weight 3 and 4 by Krawtchouk eigenvalue mass.
    For n=4, K*=5: allocates 1+12+54+54+16 = 137 operators.
    """
    from math import comb

    # Weight-class allocation for n=4, K*=5
    if n != 4 or K != 5:
        raise NotImplementedError("Only n=4, K*=5 implemented")

    allocation = {0: 1, 1: 12, 2: 54, 3: 54, 4: 16}
    rng = np.random.default_rng(0)  # deterministic selection

    operators = set()
    for w in range(n + 1):
        # All Paulis at this weight
        weight_ops = [p for p in all_paulis(n) if pauli_weight(p) == w]
        m_w = allocation[w]
        if m_w >= len(weight_ops):
            # Saturate: include all
            for p in weight_ops:
                operators.add(tuple(p))
        else:
            # Partial: select m_w deterministically
            indices = rng.choice(len(weight_ops), size=m_w, replace=False)
            for idx in indices:
                operators.add(tuple(weight_ops[idx]))

    return operators


def random_operator_indices(n, M, rng=None):
    """Return set of M random non-identity Pauli index tuples.

    Always includes the identity (0,0,...,0).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    identity = tuple([0] * n)
    all_nonid = [p for p in all_paulis(n) if p != identity]
    chosen = rng.choice(len(all_nonid), size=M - 1, replace=False)
    result = {identity}
    for idx in chosen:
        result.add(tuple(all_nonid[idx]))
    return result
