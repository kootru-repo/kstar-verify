#!/usr/bin/env python3
"""
IBM Floquet DTC K* Certification Experiment
============================================

Certifies 4-qubit Floquet discrete time-crystal (DTC) subsystems using
K*-selected Pauli measurements, demonstrating coherence structure
inaccessible to standard Z-only magnetization protocols.

Reference: "Realization of two-dimensional discrete time crystals with
anisotropic Heisenberg coupling" (Nature Comms, 2026) -- BASQ/IBM team.

Three-arm comparison:
  1. K* (137 Pauli operators) -> full density matrix reconstruction
  2. Z-only (BASQ protocol)  -> magnetization + ZZ correlators only
  3. Random (137 Paulis)      -> baseline

Usage:
  python ibm_floquet_dtc_test.py API_TOKEN [options]

Options:
  --backend NAME       IBM backend (default: least-busy Heron)
  --qubits Q0,Q1,Q2,Q3  Physical qubit indices
  --phi FLOAT          Floquet kick angle (default: pi/2, DTC phase)
  --eps FLOAT          Heisenberg anisotropy (default: 0.1)
  --J FLOAT            Base coupling strength (default: 1.0, matches BASQ)
  --disorder FLOAT     Disorder width W (default: 0.5; J_ij = J*(1+delta), delta~U[-W/2,W/2])
  --disorder-seed INT  RNG seed for disorder realization (default: 0)
  --cycles INT         Number of Floquet cycles (default: 6)
  --sweep              Multi-cycle sweep T=2,4,6,8 (main figure data)
  --n-shots INT        Shots per circuit (default: 1000)
  --n-seeds INT        Number of random seeds for baseline (default: 5)
  --init-state STATE   Initial state: neel, polarized, random (default: neel)
  --seed INT           Random seed for random Paulis (default: 42)
  --recover JOB_IDS    Recover from previous job IDs (comma-separated)
  --sim-only           Run noiseless simulation only (no IBM submission)
"""

import sys
import json
import time
import argparse
import numpy as np
from collections import defaultdict
from itertools import product as cart_product

from core import (
    select_kstar_paulis, select_random_paulis,
    state_fidelity, PAULI_LABELS, all_pauli_operators,
    I2, X, Y, Z,
    depolarizing_channel, readout_error, shot_noise,
)
from robust_mle import reconstruct_robust_mle

# Default seeds for multi-seed random baseline (5 seeds for error bars)
DEFAULT_RANDOM_SEEDS = [42, 137, 256, 512, 777]


def _make_seed_list(base_seed, n_seeds):
    """Generate n_seeds deterministic seeds from a base seed."""
    if n_seeds == 1:
        return [base_seed]
    rng = np.random.RandomState(base_seed)
    return sorted(rng.randint(1, 10000, size=n_seeds).tolist())


def _generate_bond_couplings(n_qubits, J_base=1.0, disorder_W=0.5,
                             disorder_seed=0):
    """Generate per-bond coupling strengths with disorder.

    J_ij = J_base * (1 + delta_ij) where delta_ij ~ U[-W/2, W/2].
    Matches BASQ convention (Switzer et al. Eq. 4).

    Returns:
        dict: {(i, j): J_ij} for each nearest-neighbor bond in linear chain
    """
    even_pairs = [(i, i + 1) for i in range(0, n_qubits - 1, 2)]
    odd_pairs = [(i, i + 1) for i in range(1, n_qubits - 1, 2)]
    all_pairs = even_pairs + odd_pairs

    if disorder_W == 0:
        return {pair: J_base for pair in all_pairs}

    rng = np.random.RandomState(disorder_seed)
    couplings = {}
    for pair in all_pairs:
        delta = rng.uniform(-disorder_W / 2, disorder_W / 2)
        couplings[pair] = J_base * (1 + delta)
    return couplings


# ---------------------------------------------------------------------------
# Basis grouping (qubitwise-commuting cover for circuit reduction)
# ---------------------------------------------------------------------------

def _is_compatible(pauli_string, basis_string):
    """Check if Pauli operator can be extracted from a measurement basis.
    Compatible iff for each qubit: P_i == 'I' or P_i == basis_i."""
    return all(p == 'I' or p == b for p, b in zip(pauli_string, basis_string))


def greedy_basis_cover(labels, n_qubits):
    """Find a small set of tensor-product bases covering all operators.
    Returns list of (basis_string, [covered_labels])."""
    all_bases = [''.join(b) for b in cart_product('XYZ', repeat=n_qubits)]
    uncovered = set(range(len(labels)))
    result = []
    while uncovered:
        best_basis = None
        best_cover = []
        for basis in all_bases:
            cover = [idx for idx in uncovered
                     if _is_compatible(labels[idx], basis)]
            if len(cover) > len(best_cover):
                best_cover = cover
                best_basis = basis
        if not best_basis or not best_cover:
            break
        result.append((best_basis, [labels[idx] for idx in best_cover]))
        uncovered -= set(best_cover)
    return result


def build_grouped_circuits(n_qubits, phi, eps, J, n_cycles, init_state,
                           label_sets, bond_couplings=None):
    """Build one circuit per basis for each label set (arm).

    Args:
        label_sets: list of (arm_name, pauli_labels) pairs
        bond_couplings: dict {(i,j): J_ij} for per-bond disordered couplings

    Returns:
        all_circuits: list of QuantumCircuits
        circuit_meta: list of (arm_name, basis_string, [covered_labels])
    """
    all_circuits = []
    circuit_meta = []
    for arm_name, labels in label_sets:
        cover = greedy_basis_cover(labels, n_qubits)
        for basis, covered in cover:
            qc = floquet_circuit_with_measurement(
                n_qubits, phi, eps, J, n_cycles, basis, init_state,
                bond_couplings=bond_couplings)
            qc.name = f"{arm_name}_{basis}"
            all_circuits.append(qc)
            circuit_meta.append((arm_name, basis, covered))
    return all_circuits, circuit_meta


def extract_grouped_expectations(counts_list, circuit_meta, n_qubits):
    """Extract per-operator expectations from grouped measurement results.

    Args:
        counts_list: list of count dicts (one per circuit/basis)
        circuit_meta: list of (arm_name, basis_string, [covered_labels])

    Returns:
        dict: {arm_name: {label: expectation_value}}
    """
    arm_exps = defaultdict(dict)
    for counts, (arm_name, basis, covered_labels) in zip(counts_list, circuit_meta):
        for label in covered_labels:
            arm_exps[arm_name][label] = counts_to_expectation(
                counts, label, n_qubits)
    return dict(arm_exps)


# ---------------------------------------------------------------------------
# Floquet circuit construction
# ---------------------------------------------------------------------------

def build_floquet_circuit_bare(n_qubits, phi, eps, J, n_cycles, init_state='neel',
                               bond_couplings=None):
    """
    Build Floquet DTC circuit WITHOUT measurement gates.

    Floquet unitary per cycle (matches BASQ convention, Switzer et al. Eq. 1):
        U_F = U_xxz_odd * U_xxz_even * U_kick
    Applied to |psi>: kick first, then even-bond coupling, then odd-bond.

    U_kick = prod_i exp(-i * phi * X_i)     [transverse field kick]
    U_xxz  = prod_<i,j> exp[-i*J_ij*(eps*XX + eps*YY + ZZ)]  [anisotropic Heisenberg]

    Coupling layers for linear chain: even bonds (0,1),(2,3) then odd bond (1,2).
    For n=4 this gives 3 nearest-neighbor pairs in 2 parallel layers.

    Args:
        n_qubits: number of qubits (4)
        phi: kick angle (pi/2 = DTC phase)
        eps: XY anisotropy parameter (0 = pure Ising, 1 = isotropic Heisenberg)
        J: overall coupling strength (used if bond_couplings is None)
        n_cycles: number of Floquet periods
        init_state: 'neel' (|0101>), 'polarized' (|0000>), 'afm' (|1010>)
        bond_couplings: dict {(i,j): J_ij} for per-bond disordered couplings

    Returns:
        QuantumCircuit (no measurements)
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("qiskit required for circuit construction")

    qc = QuantumCircuit(n_qubits, n_qubits)

    # --- Initial state preparation ---
    if init_state == 'neel':
        # |0101> -- Neel order (alternating, starting with 0)
        for i in range(1, n_qubits, 2):
            qc.x(i)
    elif init_state == 'afm':
        # |1010> -- anti-ferromagnetic (alternating, starting with 1)
        for i in range(0, n_qubits, 2):
            qc.x(i)
    elif init_state == 'polarized':
        pass  # |0000> is default
    else:
        raise ValueError(f"Unknown init_state: {init_state}")

    qc.barrier()

    # --- Floquet cycles ---
    # Coupling pairs for linear chain
    even_pairs = [(i, i + 1) for i in range(0, n_qubits - 1, 2)]  # (0,1), (2,3)
    odd_pairs = [(i, i + 1) for i in range(1, n_qubits - 1, 2)]   # (1,2)

    for cycle in range(n_cycles):
        # Layer 1: Transverse field kick  U_kick = prod_i Rx(2*phi)
        for i in range(n_qubits):
            qc.rx(2 * phi, i)

        # Layer 2: XXZ coupling (even bonds -- parallel)
        for (i, j) in even_pairs:
            j_ij = bond_couplings[(i, j)] if bond_couplings else J
            _append_xxz_gate(qc, i, j, j_ij, eps)

        # Layer 3: XXZ coupling (odd bonds -- parallel)
        for (i, j) in odd_pairs:
            j_ij = bond_couplings[(i, j)] if bond_couplings else J
            _append_xxz_gate(qc, i, j, j_ij, eps)

        if cycle < n_cycles - 1:
            qc.barrier()

    return qc


def _append_xxz_gate(qc, i, j, J, eps):
    """
    Append exp[-i*J*(eps*XX + eps*YY + ZZ)] to circuit on qubits (i,j).

    Decomposition using native Qiskit gates:
        Rxx(2*J*eps) * Ryy(2*J*eps) * Rzz(2*J)

    The Rxx, Ryy, Rzz gates in Qiskit implement:
        Rxx(theta) = exp(-i * theta/2 * XX)
        Ryy(theta) = exp(-i * theta/2 * YY)
        Rzz(theta) = exp(-i * theta/2 * ZZ)

    So we need theta_xx = theta_yy = 2*J*eps, theta_zz = 2*J.
    """
    qc.rxx(2 * J * eps, i, j)
    qc.ryy(2 * J * eps, i, j)
    qc.rzz(2 * J, i, j)


def floquet_circuit_with_measurement(n_qubits, phi, eps, J, n_cycles,
                                     pauli_string, init_state='neel',
                                     bond_couplings=None):
    """
    Build complete Floquet circuit with Pauli measurement basis rotation.

    Args:
        pauli_string: e.g. "XYZZ" -- which Pauli to measure on each qubit
        (other args as in build_floquet_circuit_bare)

    Returns:
        QuantumCircuit ready for submission
    """
    from qiskit import QuantumCircuit

    qc = build_floquet_circuit_bare(n_qubits, phi, eps, J, n_cycles, init_state,
                                    bond_couplings=bond_couplings)
    qc.barrier()

    # Basis rotation for Pauli measurement
    for i, basis in enumerate(pauli_string):
        if basis == 'X':
            qc.h(i)
        elif basis == 'Y':
            qc.sdg(i)
            qc.h(i)
        # Z and I: measure in computational basis (no rotation)

    qc.measure(range(n_qubits), range(n_qubits))
    return qc


# ---------------------------------------------------------------------------
# Expectation value extraction
# ---------------------------------------------------------------------------

def counts_to_expectation(counts, pauli_string, n_qubits):
    """
    Convert measurement counts to Pauli expectation value.
    Identical to ibm_hardware_test.py for consistency.
    """
    total = sum(counts.values())
    active_qubits = [i for i, p in enumerate(pauli_string) if p != 'I']

    if not active_qubits:
        return 1.0  # <I> = 1

    expectation = 0.0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]  # Qiskit little-endian
        parity = sum(int(bits[q]) for q in active_qubits) % 2
        expectation += ((-1) ** parity) * count

    return expectation / total


# ---------------------------------------------------------------------------
# Z-only order parameters (BASQ protocol)
# ---------------------------------------------------------------------------

def compute_zonly_order_params(z_expectations, zz_expectations, n_qubits,
                               init_signs, n_cycles):
    """
    Compute DTC and MBL order parameters from Z-only measurements.

    This replicates the BASQ protocol:
        Delta_DTC = (1/N) * sum_i (-1)^T * s_i * <Z_i>
        Delta_MBL = (1/N) * |sum_i s_i * <Z_i>|
        chi       = (1/M) * sum_<i,j> <Z_i Z_j>^2

    where s_i = +1/-1 for initial spin up/down, T = n_cycles.

    Args:
        z_expectations: dict {i: <Z_i>} for each qubit
        zz_expectations: dict {(i,j): <Z_i Z_j>} for nearest-neighbor pairs
        init_signs: list of +1/-1 for initial state (e.g. [+1,-1,+1,-1] for Neel)
        n_cycles: number of Floquet cycles (determines (-1)^T phase)

    Returns:
        dict with Delta_DTC, Delta_MBL, chi_spatial
    """
    N = n_qubits
    phase = (-1) ** n_cycles  # Period-doubling signature

    # Stroboscopic magnetization
    mag = sum(init_signs[i] * z_expectations[i] for i in range(N)) / N

    delta_dtc = phase * mag   # Should be positive in DTC phase
    delta_mbl = abs(mag)      # Localization persistence

    # Spatial correlator
    M = len(zz_expectations)
    if M > 0:
        chi = sum(v ** 2 for v in zz_expectations.values()) / M
    else:
        chi = 0.0

    return {
        'delta_dtc': float(delta_dtc),
        'delta_mbl': float(delta_mbl),
        'chi_spatial': float(chi),
        'stroboscopic_mag': float(mag),
    }


def extract_zonly_from_full(all_expectations, all_labels, n_qubits):
    """
    Extract Z-only observables from full K* expectation data.

    Returns:
        z_expectations: {i: <Z_i>}
        zz_expectations: {(i,j): <Z_i Z_j>}
    """
    label_to_exp = dict(zip(all_labels, all_expectations))

    # Single-qubit Z
    z_exp = {}
    for i in range(n_qubits):
        label = ['I'] * n_qubits
        label[i] = 'Z'
        label_str = ''.join(label)
        z_exp[i] = label_to_exp.get(label_str, 0.0)

    # Nearest-neighbor ZZ
    zz_exp = {}
    for i in range(n_qubits - 1):
        label = ['I'] * n_qubits
        label[i] = 'Z'
        label[i + 1] = 'Z'
        label_str = ''.join(label)
        zz_exp[(i, i + 1)] = label_to_exp.get(label_str, 0.0)

    return z_exp, zz_exp


# ---------------------------------------------------------------------------
# Coherence analysis (K*-exclusive -- invisible to Z-only)
# ---------------------------------------------------------------------------

def analyze_coherences(rho, n_qubits):
    """
    Extract coherence metrics from reconstructed density matrix that are
    invisible to Z-only measurements.

    Returns:
        dict with off-diagonal coherence metrics:
        - l1_coherence: sum of |rho_ij| for i != j (l1 norm of coherence)
        - weight_class_expectations: {w: mean |<P>| for weight-w Paulis by type}
        - x_correlators: <X_i X_j> values (invisible to Z-only)
        - y_correlators: <Y_i Y_j> values (invisible to Z-only)
        - mixed_correlators: <X_i Y_j>, <Y_i X_j> etc. (invisible to Z-only)
        - entanglement_witness: negative eigenvalue of partial transpose (if any)
    """
    dim = 2 ** n_qubits

    # l1 coherence (computational basis)
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))

    # Pauli expectations by type
    ops, labels, _ = all_pauli_operators(n_qubits)
    pauli_exps = {}
    for op, label in zip(ops, labels):
        pauli_exps[label] = np.trace(op @ rho).real

    # Weight-class breakdown
    weight_means = defaultdict(lambda: defaultdict(list))
    for label, exp in pauli_exps.items():
        if label == 'I' * n_qubits:
            continue
        weight = sum(1 for c in label if c != 'I')
        # Classify by Pauli type composition
        has_x = 'X' in label
        has_y = 'Y' in label
        has_z = 'Z' in label
        if has_z and not has_x and not has_y:
            ptype = 'Z-only'
        elif (has_x or has_y) and not has_z:
            ptype = 'XY-only'
        else:
            ptype = 'mixed'
        weight_means[weight][ptype].append(abs(exp))

    weight_summary = {}
    for w in sorted(weight_means.keys()):
        weight_summary[w] = {}
        for ptype, vals in weight_means[w].items():
            weight_summary[w][ptype] = {
                'mean_abs': float(np.mean(vals)),
                'max_abs': float(np.max(vals)),
                'count': len(vals),
            }

    # XX, YY correlators (nearest-neighbor)
    xx_corr = {}
    yy_corr = {}
    xy_corr = {}
    for i in range(n_qubits - 1):
        j = i + 1
        # XX
        label = ['I'] * n_qubits
        label[i] = 'X'
        label[j] = 'X'
        xx_corr[(i, j)] = pauli_exps.get(''.join(label), 0.0)
        # YY
        label = ['I'] * n_qubits
        label[i] = 'Y'
        label[j] = 'Y'
        yy_corr[(i, j)] = pauli_exps.get(''.join(label), 0.0)
        # XY
        label = ['I'] * n_qubits
        label[i] = 'X'
        label[j] = 'Y'
        xy_corr[(i, j)] = pauli_exps.get(''.join(label), 0.0)

    # Partial transpose entanglement witness (bipartition 0,1 | 2,3)
    pt_neg = _partial_transpose_negativity(rho, n_qubits)

    return {
        'l1_coherence': float(l1),
        'weight_class_summary': weight_summary,
        'xx_correlators': {f"{k[0]},{k[1]}": float(v) for k, v in xx_corr.items()},
        'yy_correlators': {f"{k[0]},{k[1]}": float(v) for k, v in yy_corr.items()},
        'xy_correlators': {f"{k[0]},{k[1]}": float(v) for k, v in xy_corr.items()},
        'partial_transpose_negativity': float(pt_neg),
    }


def _partial_transpose_negativity(rho, n_qubits):
    """
    Compute negativity of the partial transpose across the middle bipartition.
    Negativity > 0 implies entanglement (PPT criterion).
    """
    dim = 2 ** n_qubits
    n_a = n_qubits // 2
    n_b = n_qubits - n_a
    d_a = 2 ** n_a
    d_b = 2 ** n_b

    # Reshape into (d_a, d_b, d_a, d_b), partial transpose on subsystem A
    rho_reshaped = rho.reshape(d_a, d_b, d_a, d_b)
    rho_pt = rho_reshaped.transpose(2, 1, 0, 3).reshape(dim, dim)

    eigvals = np.linalg.eigvalsh(rho_pt)
    negativity = sum(abs(e) for e in eigvals if e < -1e-10)
    return negativity


# ---------------------------------------------------------------------------
# K*-accessible coherence analysis (C5: from expectations, not rho)
# ---------------------------------------------------------------------------

def analyze_kstar_accessible_coherences(expectations, labels, n_qubits):
    """
    Extract coherence metrics directly from K* expectation values,
    WITHOUT reconstructing the density matrix. These are the observables
    that K* actually measures -- model-independent, no reference state needed.

    Returns:
        dict with metrics computable directly from expectations
    """
    label_to_exp = dict(zip(labels, expectations))

    # Classify expectations by type
    z_only_exps = []
    xy_exps = []
    mixed_exps = []

    for label, exp in label_to_exp.items():
        if label == 'I' * n_qubits:
            continue
        has_x = 'X' in label
        has_y = 'Y' in label
        has_z = 'Z' in label
        if has_z and not has_x and not has_y:
            z_only_exps.append(abs(exp))
        elif (has_x or has_y) and not has_z:
            xy_exps.append(abs(exp))
        else:
            mixed_exps.append(abs(exp))

    # Nearest-neighbor correlators (directly measured, no reconstruction)
    xx_corr = {}
    yy_corr = {}
    xy_corr = {}
    for i in range(n_qubits - 1):
        j = i + 1
        l_xx = ['I'] * n_qubits; l_xx[i] = 'X'; l_xx[j] = 'X'
        l_yy = ['I'] * n_qubits; l_yy[i] = 'Y'; l_yy[j] = 'Y'
        l_xy = ['I'] * n_qubits; l_xy[i] = 'X'; l_xy[j] = 'Y'
        xx_corr[(i, j)] = label_to_exp.get(''.join(l_xx), 0.0)
        yy_corr[(i, j)] = label_to_exp.get(''.join(l_yy), 0.0)
        xy_corr[(i, j)] = label_to_exp.get(''.join(l_xy), 0.0)

    # Signal-to-noise: any nonzero XY expectation is a coherence detection
    xy_signal = sum(xy_exps) if xy_exps else 0.0
    total_signal = sum(z_only_exps) + xy_signal + sum(mixed_exps)
    xy_fraction = xy_signal / total_signal if total_signal > 0 else 0.0

    return {
        'n_z_only': len(z_only_exps),
        'n_xy': len(xy_exps),
        'n_mixed': len(mixed_exps),
        'mean_abs_z_only': float(np.mean(z_only_exps)) if z_only_exps else 0.0,
        'mean_abs_xy': float(np.mean(xy_exps)) if xy_exps else 0.0,
        'mean_abs_mixed': float(np.mean(mixed_exps)) if mixed_exps else 0.0,
        'xy_signal_fraction': float(xy_fraction),
        'xx_correlators': {f"{k[0]},{k[1]}": float(v) for k, v in xx_corr.items()},
        'yy_correlators': {f"{k[0]},{k[1]}": float(v) for k, v in yy_corr.items()},
        'xy_correlators': {f"{k[0]},{k[1]}": float(v) for k, v in xy_corr.items()},
        'any_xy_nonzero': any(abs(e) > 0.05 for e in xy_exps),
    }


# ---------------------------------------------------------------------------
# Noisy simulation (C4: depolarizing + readout error model)
# ---------------------------------------------------------------------------

def simulate_noisy_expectations(ideal_exps, labels, n_shots,
                                depol_rate=0.02, readout_rate=0.02):
    """
    Apply a hardware-realistic noise model to ideal expectations.

    Pipeline: ideal -> depolarizing attenuation -> readout error -> shot noise

    This is more conservative than noiseless + shot noise alone.

    Args:
        ideal_exps: dict {label: <P>_ideal}
        labels: list of Pauli labels to extract
        n_shots: shot count for binomial noise
        depol_rate: per-qubit depolarizing probability
        readout_rate: per-qubit readout error probability

    Returns:
        list of noisy expectation values in label order
    """
    noisy = []
    for label in labels:
        exp = ideal_exps.get(label, 0.0)
        # Depolarizing: weight-w Pauli attenuated by (1-2p)^w
        w = sum(1 for c in label if c != 'I')
        exp *= (1 - 2 * depol_rate) ** w
        # Readout error: independent per-qubit, attenuation (1-2*readout_rate)^w
        exp *= (1 - 2 * readout_rate) ** w
        # Shot noise
        exp = shot_noise(exp, n_shots)
        noisy.append(exp)
    return noisy


# ---------------------------------------------------------------------------
# Qubit ordering cross-validation (C1)
# ---------------------------------------------------------------------------

def verify_qubit_ordering(verbose=True):
    """
    Cross-validate that the numpy simulation and Qiskit circuit convention
    produce consistent Pauli expectations. Tests on Neel state at T=0
    (no evolution) where expectations are analytically known.

    Neel |0101> expectations:
        <Z_0> = +1, <Z_1> = -1, <Z_2> = +1, <Z_3> = -1
        <Z_0 Z_1> = -1, <Z_1 Z_2> = -1, <Z_2 Z_3> = -1
        <ZZZZ> = +1
        All X, Y single-qubit: 0

    Returns True if all checks pass.
    """
    n = 4
    # T=0: just the initial state (use T=0 which means no Floquet cycles)
    # simulate_floquet_noiseless with n_cycles=0 should give the bare Neel state
    rho = simulate_floquet_noiseless(n, np.pi / 2, 0.1, 0.8, 0, 'neel')

    ops, labels, _ = all_pauli_operators(n)
    exps = {l: np.trace(o @ rho).real for o, l in zip(ops, labels)}

    checks = []

    # Single-qubit Z: alternating +1, -1
    expected_z = {0: +1, 1: -1, 2: +1, 3: -1}
    for i, ez in expected_z.items():
        l = ['I'] * n; l[i] = 'Z'
        label = ''.join(l)
        val = exps[label]
        ok = abs(val - ez) < 1e-10
        checks.append(ok)
        if verbose and not ok:
            print(f"  FAIL: <{label}> = {val:.6f}, expected {ez}")

    # ZZ nearest-neighbor: all -1 for Neel
    for i in range(n - 1):
        l = ['I'] * n; l[i] = 'Z'; l[i + 1] = 'Z'
        label = ''.join(l)
        val = exps[label]
        ok = abs(val - (-1.0)) < 1e-10
        checks.append(ok)
        if verbose and not ok:
            print(f"  FAIL: <{label}> = {val:.6f}, expected -1")

    # ZZZZ = +1 (product of all Z eigenvalues: +1*-1*+1*-1 = +1)
    val = exps['ZZZZ']
    ok = abs(val - 1.0) < 1e-10
    checks.append(ok)
    if verbose and not ok:
        print(f"  FAIL: <ZZZZ> = {val:.6f}, expected +1")

    # All single-qubit X and Y should be 0
    for i in range(n):
        for basis in ['X', 'Y']:
            l = ['I'] * n; l[i] = basis
            label = ''.join(l)
            val = exps[label]
            ok = abs(val) < 1e-10
            checks.append(ok)
            if verbose and not ok:
                print(f"  FAIL: <{label}> = {val:.6f}, expected 0")

    passed = all(checks)
    if verbose:
        print(f"  Qubit ordering verification: {'PASS' if passed else 'FAIL'} "
              f"({sum(checks)}/{len(checks)} checks)")
    return passed


# ---------------------------------------------------------------------------
# Multi-seed random reconstruction helper
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Noiseless simulation (for validation and sim-only mode)
# ---------------------------------------------------------------------------

def _init_psi(n_qubits, init_state):
    """Prepare initial statevector. Shared by simulation functions."""
    dim = 2 ** n_qubits
    psi = np.zeros(dim, dtype=complex)
    if init_state == 'neel':
        idx = 0
        for i in range(n_qubits):
            if i % 2 == 1:
                idx |= (1 << (n_qubits - 1 - i))
        psi[idx] = 1.0
    elif init_state == 'afm':
        idx = 0
        for i in range(n_qubits):
            if i % 2 == 0:
                idx |= (1 << (n_qubits - 1 - i))
        psi[idx] = 1.0
    elif init_state == 'polarized':
        psi[0] = 1.0
    return psi


def _build_floquet_unitary(n_qubits, phi, eps, J, bond_couplings=None):
    """Build the full Floquet unitary. Shared by simulation functions.

    Ordering matches BASQ (Switzer et al. Eq. 1): kick first, then
    even-bond coupling, then odd-bond coupling.
    """
    dim = 2 ** n_qubits

    rx = np.array([
        [np.cos(phi), -1j * np.sin(phi)],
        [-1j * np.sin(phi), np.cos(phi)]
    ], dtype=complex)
    U_kick = rx
    for _ in range(n_qubits - 1):
        U_kick = np.kron(U_kick, rx)

    def xxz_unitary_2q(J_val, eps_val):
        H_xxz = J_val * (eps_val * np.kron(X, X) +
                         eps_val * np.kron(Y, Y) +
                         np.kron(Z, Z))
        from scipy.linalg import expm
        return expm(-1j * H_xxz)

    def embed_2q(U_2q_local, q0, q1, n):
        dim_full = 2 ** n
        U_full = np.eye(dim_full, dtype=complex)
        for row in range(dim_full):
            for col in range(dim_full):
                match = True
                for q in range(n):
                    if q == q0 or q == q1:
                        continue
                    if ((row >> (n - 1 - q)) & 1) != ((col >> (n - 1 - q)) & 1):
                        match = False
                        break
                if not match:
                    continue
                r0 = (row >> (n - 1 - q0)) & 1
                r1 = (row >> (n - 1 - q1)) & 1
                c0 = (col >> (n - 1 - q0)) & 1
                c1 = (col >> (n - 1 - q1)) & 1
                U_full[row, col] = U_2q_local[r0 * 2 + r1, c0 * 2 + c1]
        return U_full

    even_pairs = [(i, i + 1) for i in range(0, n_qubits - 1, 2)]
    odd_pairs = [(i, i + 1) for i in range(1, n_qubits - 1, 2)]

    # Per-bond unitaries (supports disordered J_ij)
    U_even_list = []
    for (q0, q1) in even_pairs:
        j_ij = bond_couplings[(q0, q1)] if bond_couplings else J
        U_even_list.append(embed_2q(xxz_unitary_2q(j_ij, eps), q0, q1, n_qubits))

    U_odd_list = []
    for (q0, q1) in odd_pairs:
        j_ij = bond_couplings[(q0, q1)] if bond_couplings else J
        U_odd_list.append(embed_2q(xxz_unitary_2q(j_ij, eps), q0, q1, n_qubits))

    U_even = np.eye(dim, dtype=complex)
    for U_e in U_even_list:
        U_even = U_e @ U_even

    U_odd = np.eye(dim, dtype=complex)
    for U_o in U_odd_list:
        U_odd = U_o @ U_odd

    return U_odd @ U_even @ U_kick


def simulate_floquet_noiseless(n_qubits, phi, eps, J, n_cycles, init_state='neel',
                               bond_couplings=None):
    """
    Exact statevector simulation of the Floquet circuit.

    Args:
        bond_couplings: dict {(i,j): J_ij} for per-bond disordered couplings.
                        If None, uses uniform J on all bonds.

    Returns:
        rho_ideal: density matrix after n_cycles Floquet periods
    """
    psi = _init_psi(n_qubits, init_state)
    U_floquet = _build_floquet_unitary(n_qubits, phi, eps, J, bond_couplings)

    for _ in range(n_cycles):
        psi = U_floquet @ psi

    return np.outer(psi, psi.conj())


# ---------------------------------------------------------------------------
# Stroboscopic time trace (P4: per-cycle Delta(t), BASQ Eq. 5)
# ---------------------------------------------------------------------------

def compute_stroboscopic_trace(n_qubits, phi, eps, J, n_cycles, init_state='neel',
                               bond_couplings=None):
    """Compute Delta(t) for t=0,1,...,n_cycles (BASQ Eq. 5).

    Delta(t) = (1/N) sum_i s_i <Z_i>(t)
    where s_i are the initial spin signs.

    Returns:
        list of dicts: [{'t': t, 'delta': Delta(t), 'delta_dtc': (-1)^t * Delta(t)}]
    """
    psi = _init_psi(n_qubits, init_state)
    U_floquet = _build_floquet_unitary(n_qubits, phi, eps, J, bond_couplings)
    init_signs = _init_signs_for_state(init_state, n_qubits)
    dim = 2 ** n_qubits

    # Build diagonal Z operators for each qubit
    Z_diags = []
    for i in range(n_qubits):
        diag = np.array([(-1) ** ((idx >> (n_qubits - 1 - i)) & 1)
                         for idx in range(dim)], dtype=float)
        Z_diags.append(diag)

    trace = []
    for t in range(n_cycles + 1):
        if t > 0:
            psi = U_floquet @ psi
        probs = np.abs(psi) ** 2
        delta_t = sum(init_signs[i] * np.dot(Z_diags[i], probs)
                      for i in range(n_qubits)) / n_qubits
        trace.append({
            't': t,
            'delta': float(delta_t),
            'delta_dtc': float((-1) ** t * delta_t),
        })
    return trace


# ---------------------------------------------------------------------------
# Hamming distance distribution (P2: from rho diagonal)
# ---------------------------------------------------------------------------

def compute_hamming_distribution(rho, init_state, n_qubits):
    """Compute Hamming distance distribution from the density matrix diagonal.

    Hamming distance d = number of bit flips from the initial state.
    Matches BASQ Fig. 2c,d observable.

    Returns:
        dict with 'distribution' (list of {distance, probability}),
        'mean_distance', 'variance'
    """
    dim = 2 ** n_qubits
    populations = np.real(np.diag(rho))

    # Determine initial bitstring
    if init_state == 'neel':
        init_bits = [0 if i % 2 == 0 else 1 for i in range(n_qubits)]
    elif init_state == 'afm':
        init_bits = [1 if i % 2 == 0 else 0 for i in range(n_qubits)]
    elif init_state == 'polarized':
        init_bits = [0] * n_qubits
    else:
        init_bits = [0] * n_qubits

    # Compute Hamming distance for each basis state
    dist = defaultdict(float)
    for idx in range(dim):
        bits = [(idx >> (n_qubits - 1 - i)) & 1 for i in range(n_qubits)]
        hamming_d = sum(b != ib for b, ib in zip(bits, init_bits))
        dist[hamming_d] += populations[idx]

    distribution = [{'distance': d, 'probability': float(dist[d])}
                    for d in range(n_qubits + 1)]
    mean_d = sum(d * dist[d] for d in range(n_qubits + 1))
    var_d = sum(d ** 2 * dist[d] for d in range(n_qubits + 1)) - mean_d ** 2

    return {
        'distribution': distribution,
        'mean_distance': float(mean_d),
        'variance': float(var_d),
    }


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _extract_all_expectations(all_results, circuit_labels, n_qubits):
    """Extract expectation values from IBM job results (V2 API)."""
    all_expectations = []
    for batch_start, result in all_results:
        for j in range(len(result)):
            pub_result = result[j]
            counts = None
            try:
                counts = pub_result.data.meas.get_counts()
            except AttributeError:
                data_bin = pub_result.data
                for attr_name in dir(data_bin):
                    if attr_name.startswith('_'):
                        continue
                    attr = getattr(data_bin, attr_name)
                    if hasattr(attr, 'get_counts'):
                        counts = attr.get_counts()
                        break
            if counts is None:
                all_expectations.append(0.0)
                continue
            idx = batch_start + j
            label = circuit_labels[idx]
            exp = counts_to_expectation(counts, label, n_qubits)
            all_expectations.append(exp)
    return all_expectations


def _init_signs_for_state(init_state, n_qubits):
    """Return initial spin signs for order parameter computation."""
    if init_state == 'neel':
        return [+1 if i % 2 == 0 else -1 for i in range(n_qubits)]
    elif init_state == 'afm':
        return [-1 if i % 2 == 0 else +1 for i in range(n_qubits)]
    return [+1] * n_qubits


# ---------------------------------------------------------------------------
# IBM hardware execution
# ---------------------------------------------------------------------------

def run_floquet_dtc_test(api_token, n_qubits=4, n_shots=1000,
                         backend_name=None, initial_qubits=None,
                         phi=np.pi / 2, eps=0.1, J=1.0, n_cycles=6,
                         init_state='neel', seed=42, n_seeds=5,
                         full_tomo=False, grouped=False,
                         disorder_W=0.5, disorder_seed=0,
                         instance=None):
    """
    Multi-arm Floquet DTC certification on IBM hardware.

    Arms:
        1. K* (137 Paulis) -> robust MLE -> full density matrix
        2. Z-only (BASQ protocol) -> order parameters only (extracted from K*)
        3. Random x n_seeds (137 Paulis each) -> robust MLE -> mean +- std
        4. (optional) Full tomo (255 Paulis) -> gold-standard reference

    Total circuits: 137 + n_seeds*137 [+ 255 if full_tomo].
    """
    from qiskit import QuantumCircuit
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    from qiskit.transpiler import preset_passmanagers

    seeds = _make_seed_list(seed, n_seeds)
    bond_couplings = _generate_bond_couplings(n_qubits, J, disorder_W, disorder_seed)

    print("=" * 70)
    print("  FLOQUET DTC K* CERTIFICATION EXPERIMENT")
    print("=" * 70)
    print(f"  Qubits:      {n_qubits}")
    print(f"  Floquet:     phi={phi:.4f}, eps={eps}, J={J}, T={n_cycles}")
    print(f"  Disorder:    W={disorder_W}, seed={disorder_seed}")
    bc_str = ', '.join(f"({i},{j}):{v:.3f}" for (i, j), v in sorted(bond_couplings.items()))
    print(f"  Bond J_ij:   {bc_str}")
    print(f"  Init state:  {init_state}")
    print(f"  Shots:       {n_shots}")
    print(f"  Random seeds: {seeds}")
    if full_tomo:
        print(f"  Full tomo:   YES (gold-standard arm)")
    print()

    # --- Step 1: Select Pauli operators ---
    kstar_ops, kstar_labels, _ = select_kstar_paulis(n_qubits)
    n_kstar = len(kstar_labels)

    # Multiple random seeds
    rand_arms = []
    for s in seeds:
        r_ops, r_labels, _ = select_random_paulis(n_qubits, n_kstar, seed=s)
        rand_arms.append((s, r_ops, r_labels))

    # Full tomo (optional)
    full_tomo_labels = None
    full_tomo_ops = None
    if full_tomo:
        from itertools import product as iterproduct
        full_tomo_labels = [''.join(p) for p in iterproduct('IXYZ', repeat=n_qubits)]
        full_tomo_labels = [l for l in full_tomo_labels if l != 'I' * n_qubits]
        full_tomo_ops_all, full_tomo_labels_all, _ = all_pauli_operators(n_qubits)
        # Build ops list matching labels order
        label_to_op = dict(zip(full_tomo_labels_all, full_tomo_ops_all))
        full_tomo_ops = [label_to_op[l] for l in full_tomo_labels]

    n_rand_total = n_seeds * n_kstar
    n_full = len(full_tomo_labels) if full_tomo else 0
    total_circuits = n_kstar + n_rand_total + n_full

    print(f"  K* Paulis:     {n_kstar}")
    print(f"  Random Paulis: {n_kstar} x {n_seeds} seeds = {n_rand_total}")
    if full_tomo:
        print(f"  Full tomo:     {n_full}")
    print(f"  Total circuits: {total_circuits}")
    print()

    # --- Step 2: Build circuits ---
    print("  Building Floquet + measurement circuits...")

    if grouped:
        # Basis-grouped mode: one circuit per basis, not per operator
        label_sets = [('kstar', list(kstar_labels))]
        for s, r_ops, r_labels in rand_arms:
            label_sets.append((f'random_s{s}', list(r_labels)))
        if full_tomo:
            label_sets.append(('full_tomo', list(full_tomo_labels)))

        all_circuits, circuit_meta = build_grouped_circuits(
            n_qubits, phi, eps, J, n_cycles, init_state, label_sets,
            bond_couplings=bond_couplings)
        circuit_labels = [basis for _, basis, _ in circuit_meta]

        n_bases_kstar = sum(1 for arm, _, _ in circuit_meta if arm == 'kstar')
        n_bases_rand = sum(1 for arm, _, _ in circuit_meta if arm.startswith('random'))
        total_circuits = len(all_circuits)
        print(f"  GROUPED MODE: {n_bases_kstar} K* bases + {n_bases_rand} random bases"
              f" = {total_circuits} circuits")
        print(f"  (vs {n_kstar + n_rand_total + n_full} ungrouped)")
    else:
        all_circuits = []
        circuit_labels = []
        circuit_meta = None  # not used in ungrouped mode

        # K* arm
        for label in kstar_labels:
            qc = floquet_circuit_with_measurement(
                n_qubits, phi, eps, J, n_cycles, label, init_state,
                bond_couplings=bond_couplings)
            qc.name = f"kstar_{label}"
            all_circuits.append(qc)
            circuit_labels.append(label)

        # Random arms (per seed)
        for s, r_ops, r_labels in rand_arms:
            for label in r_labels:
                qc = floquet_circuit_with_measurement(
                    n_qubits, phi, eps, J, n_cycles, label, init_state,
                    bond_couplings=bond_couplings)
                qc.name = f"rand_s{s}_{label}"
                all_circuits.append(qc)
                circuit_labels.append(label)

        # Full tomo arm
        if full_tomo:
            for label in full_tomo_labels:
                qc = floquet_circuit_with_measurement(
                    n_qubits, phi, eps, J, n_cycles, label, init_state,
                    bond_couplings=bond_couplings)
                qc.name = f"full_{label}"
                all_circuits.append(qc)
                circuit_labels.append(label)

        total_circuits = len(all_circuits)

    depth_sample = all_circuits[0].depth()
    print(f"  Circuit depth (pre-transpile, sample): {depth_sample}")
    print()

    # --- Step 3: Connect to IBM ---
    print("  Connecting to IBM Quantum...")
    svc_kwargs = dict(channel="ibm_quantum_platform", token=api_token)
    if instance:
        svc_kwargs['instance'] = instance
    service = QiskitRuntimeService(**svc_kwargs)

    if backend_name:
        backend = service.backend(backend_name)
    else:
        backends = service.backends(
            filters=lambda b: b.num_qubits >= n_qubits and b.status().operational
        )
        backend = min(backends, key=lambda b: b.status().pending_jobs)

    print(f"  Backend:     {backend.name}")
    print(f"  Pending:     {backend.status().pending_jobs} jobs")
    print()

    # --- Step 4: Transpile ---
    print("  Transpiling...")
    pm_kwargs = {'optimization_level': 1, 'backend': backend}
    if initial_qubits:
        pm_kwargs['initial_layout'] = initial_qubits
    pm = preset_passmanagers.generate_preset_pass_manager(**pm_kwargs)
    transpiled = pm.run(all_circuits)

    transpiled_depth = transpiled[0].depth()
    print(f"  Transpiled depth (sample): {transpiled_depth}")
    print()

    # --- Step 5: Submit in batches ---
    print("  Submitting to IBM...")
    batch_size = 300 if grouped else 100
    all_results = []
    job_ids = []

    sampler = SamplerV2(mode=backend)
    for i in range(0, len(transpiled), batch_size):
        batch = transpiled[i:i + batch_size]
        batch_pubs = [(circ,) for circ in batch]
        job = sampler.run(batch_pubs)
        job_id = job.job_id()
        job_ids.append(job_id)
        print(f"    Batch {i // batch_size + 1}: {len(batch)} circuits, job={job_id}")

        t0 = time.time()
        result = job.result()
        elapsed = time.time() - t0
        print(f"    -> completed in {elapsed:.1f}s")
        all_results.append((i, result))

    # Save job IDs for recovery
    recovery_file = f"floquet_jobs_{backend.name}_{time.strftime('%Y%m%d_%H%M%S')}.json"
    with open(recovery_file, 'w') as f:
        json.dump({
            'job_ids': job_ids,
            'backend': backend.name,
            'params': {'phi': phi, 'eps': eps, 'J': J, 'n_cycles': n_cycles,
                       'init_state': init_state, 'n_qubits': n_qubits,
                       'n_shots': n_shots, 'seed': seed,
                       'disorder_W': disorder_W, 'disorder_seed': disorder_seed},
            'kstar_labels': list(kstar_labels),
            'n_seeds': n_seeds,
            'random_seeds': seeds,
            'rand_arms': [{'seed': s, 'labels': list(rl)}
                          for s, _, rl in rand_arms],
            'full_tomo': full_tomo,
        }, f, indent=2)
    print(f"  Recovery file: {recovery_file}")
    print()

    # --- Step 6: Extract expectations ---
    print("  Extracting expectation values...")

    if grouped and circuit_meta is not None:
        # Grouped mode: extract counts per circuit, then fan out per operator
        all_counts = []
        for batch_start, result in all_results:
            for j in range(len(result)):
                pub_result = result[j]
                counts = None
                try:
                    counts = pub_result.data.meas.get_counts()
                except AttributeError:
                    data_bin = pub_result.data
                    for attr_name in dir(data_bin):
                        if attr_name.startswith('_'):
                            continue
                        attr = getattr(data_bin, attr_name)
                        if hasattr(attr, 'get_counts'):
                            counts = attr.get_counts()
                            break
                all_counts.append(counts if counts else {})

        arm_exps = extract_grouped_expectations(all_counts, circuit_meta, n_qubits)
        kstar_exp_dict = arm_exps.get('kstar', {})
        kstar_expectations = [kstar_exp_dict.get(l, 0.0) for l in kstar_labels]
        rand_seed_expectations = []
        for s, _, r_labels in rand_arms:
            arm_key = f'random_s{s}'
            rd = arm_exps.get(arm_key, {})
            rand_seed_expectations.append([rd.get(l, 0.0) for l in r_labels])
        full_tomo_expectations = None
        if full_tomo:
            ftd = arm_exps.get('full_tomo', {})
            full_tomo_expectations = [ftd.get(l, 0.0) for l in full_tomo_labels]
    else:
        all_expectations = _extract_all_expectations(all_results, circuit_labels, n_qubits)
        # Slice into arms
        idx = 0
        kstar_expectations = all_expectations[idx:idx + n_kstar]; idx += n_kstar
        rand_seed_expectations = []
        for s, r_ops, r_labels in rand_arms:
            rand_seed_expectations.append(all_expectations[idx:idx + n_kstar])
            idx += n_kstar
        full_tomo_expectations = None
        if full_tomo:
            full_tomo_expectations = all_expectations[idx:idx + n_full]

    # --- Step 7: Reconstruct density matrices ---
    print("  Reconstructing density matrices (robust MLE)...")

    rho_kstar = reconstruct_robust_mle(
        kstar_expectations, kstar_ops, n_qubits, n_shots=n_shots)

    # Multi-seed random reconstruction
    rand_rhos = []
    for i, (s, r_ops, r_labels) in enumerate(rand_arms):
        rho_r = reconstruct_robust_mle(
            rand_seed_expectations[i], r_ops, n_qubits, n_shots=n_shots)
        rand_rhos.append((s, rho_r))

    # Full tomo reconstruction (gold standard)
    rho_full = None
    if full_tomo:
        rho_full = reconstruct_robust_mle(
            full_tomo_expectations, full_tomo_ops, n_qubits, n_shots=n_shots)

    # --- Step 8: Compute ideal-state fidelities ---
    print("  Running noiseless simulation for ideal-state reference...")
    rho_ideal = simulate_floquet_noiseless(n_qubits, phi, eps, J, n_cycles, init_state,
                                           bond_couplings=bond_couplings)

    f_kstar = state_fidelity(rho_ideal, rho_kstar)
    f_rands = [state_fidelity(rho_ideal, rho_r) for (s, rho_r) in rand_rhos]
    f_rand_mean = float(np.mean(f_rands))
    f_rand_std = float(np.std(f_rands))
    advantage = f_kstar - f_rand_mean

    print()
    print("  === IDEAL-STATE FIDELITY RESULTS ===")
    print(f"  F(K*)          = {f_kstar:.4f}")
    print(f"  F(rand) mean   = {f_rand_mean:.4f} +- {f_rand_std:.4f}  (n={n_seeds} seeds)")
    for s, fr in zip(seeds, f_rands):
        print(f"    seed {s:>5d}: F = {fr:.4f}")
    print(f"  dF             = {advantage:+.4f}")
    if full_tomo and rho_full is not None:
        f_full = state_fidelity(rho_ideal, rho_full)
        print(f"  F(full tomo)   = {f_full:.4f}  (gold standard, {n_full} Paulis)")
    print()

    # --- Step 9: Z-only order parameters (extracted from K* data) ---
    print("  Computing Z-only order parameters...")
    z_exp, zz_exp = extract_zonly_from_full(kstar_expectations, kstar_labels, n_qubits)
    init_signs = _init_signs_for_state(init_state, n_qubits)
    zonly_params = compute_zonly_order_params(z_exp, zz_exp, n_qubits,
                                              init_signs, n_cycles)

    print(f"  Delta_DTC = {zonly_params['delta_dtc']:.4f}")
    print(f"  Delta_MBL = {zonly_params['delta_mbl']:.4f}")
    print(f"  chi_spatial = {zonly_params['chi_spatial']:.4f}")
    print()

    # --- Step 10: Coherence analysis ---
    print("  Analyzing coherences...")

    # K*-accessible (model-independent, from expectations directly)
    kstar_direct_coh = analyze_kstar_accessible_coherences(
        kstar_expectations, kstar_labels, n_qubits)

    # Full coherence from reconstructed rho
    coherence_kstar = analyze_coherences(rho_kstar, n_qubits)
    coherence_ideal = analyze_coherences(rho_ideal, n_qubits)

    print(f"  K*-direct XY signal fraction: {kstar_direct_coh['xy_signal_fraction']:.4f}")
    print(f"  K*-direct any XY nonzero:     {kstar_direct_coh['any_xy_nonzero']}")
    print(f"  l1 coherence (ideal):  {coherence_ideal['l1_coherence']:.4f}")
    print(f"  l1 coherence (K*):     {coherence_kstar['l1_coherence']:.4f}")
    print(f"  PT negativity (ideal): {coherence_ideal['partial_transpose_negativity']:.4f}")
    print(f"  PT negativity (K*):    {coherence_kstar['partial_transpose_negativity']:.4f}")
    print()

    # --- Step 11: Self-consistency metrics (model-independent) ---
    rand_rho_mean = np.mean([rho_r for (s, rho_r) in rand_rhos], axis=0)
    trace_dist_kr = 0.5 * np.linalg.norm(rho_kstar - rand_rho_mean, 'nuc')

    # --- Step 11b: Stroboscopic trace and Hamming distribution ---
    print("  Computing stroboscopic trace and Hamming distribution...")
    strobe_trace = compute_stroboscopic_trace(
        n_qubits, phi, eps, J, n_cycles, init_state, bond_couplings)
    hamming_kstar = compute_hamming_distribution(rho_kstar, init_state, n_qubits)
    hamming_ideal = compute_hamming_distribution(rho_ideal, init_state, n_qubits)

    # --- Step 12: Save results ---
    results = {
        'experiment': 'floquet_dtc_kstar_certification',
        'backend': backend.name,
        'n_qubits': n_qubits,
        'n_shots': n_shots,
        'n_seeds': n_seeds,
        'random_seeds': seeds,
        'floquet_params': {
            'phi': float(phi), 'eps': float(eps), 'J': float(J),
            'n_cycles': n_cycles, 'init_state': init_state,
            'disorder_W': float(disorder_W), 'disorder_seed': int(disorder_seed),
            'bond_couplings': {f"{k[0]},{k[1]}": float(v)
                               for k, v in bond_couplings.items()},
        },
        'circuit_info': {
            'depth_pretranspile': depth_sample,
            'depth_transpiled': transpiled_depth,
            'n_kstar_paulis': n_kstar,
            'n_rand_paulis_per_seed': n_kstar,
            'n_rand_seeds': n_seeds,
            'n_full_tomo_paulis': n_full,
            'total_circuits': total_circuits,
        },
        'fidelity': {
            'note': 'ideal-state fidelity (vs noiseless simulation, not ground truth)',
            'f_kstar': float(f_kstar),
            'f_rand_mean': f_rand_mean,
            'f_rand_std': f_rand_std,
            'f_rand_per_seed': {str(s): float(fr) for s, fr in zip(seeds, f_rands)},
            'advantage_mean': float(advantage),
        },
        'self_consistency': {
            'trace_distance_kstar_vs_rand_mean': float(trace_dist_kr),
        },
        'zonly_order_params': zonly_params,
        'kstar_direct_coherences': {
            'xy_signal_fraction': kstar_direct_coh['xy_signal_fraction'],
            'any_xy_nonzero': kstar_direct_coh['any_xy_nonzero'],
            'xx_correlators': kstar_direct_coh['xx_correlators'],
            'yy_correlators': kstar_direct_coh['yy_correlators'],
            'xy_correlators': kstar_direct_coh['xy_correlators'],
        },
        'coherence_from_rho': {
            'kstar': {
                'l1_coherence': coherence_kstar['l1_coherence'],
                'partial_transpose_negativity': coherence_kstar['partial_transpose_negativity'],
            },
            'ideal': {
                'l1_coherence': coherence_ideal['l1_coherence'],
                'partial_transpose_negativity': coherence_ideal['partial_transpose_negativity'],
            },
        },
        'kstar_expectations': [float(e) for e in kstar_expectations],
        'kstar_labels': list(kstar_labels),
        'rand_expectations_per_seed': {
            str(s): [float(e) for e in rand_seed_expectations[i]]
            for i, (s, _, _) in enumerate(rand_arms)
        },
        'z_expectations': {str(k): float(v) for k, v in z_exp.items()},
        'zz_expectations': {f"{k[0]},{k[1]}": float(v) for k, v in zz_exp.items()},
        'stroboscopic_trace': strobe_trace,
        'hamming_distribution': {
            'kstar': hamming_kstar,
            'ideal': hamming_ideal,
        },
        'density_matrices': {
            'kstar': {'real': rho_kstar.real.tolist(), 'imag': rho_kstar.imag.tolist()},
            'ideal': {'real': rho_ideal.real.tolist(), 'imag': rho_ideal.imag.tolist()},
        },
        'job_ids': job_ids,
        'physical_qubits': initial_qubits,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    }
    if full_tomo and rho_full is not None:
        results['fidelity']['f_full_tomo'] = float(f_full)

    outfile = f"floquet_dtc_results_{backend.name}_{time.strftime('%Y%m%d_%H%M%S')}.json"
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"  Results saved: {outfile}")
    print()

    # --- Summary ---
    print("=" * 70)
    print("  SUMMARY: Floquet DTC K* Certification")
    print("=" * 70)
    print(f"  F(K*)  = {f_kstar:.4f}")
    print(f"  F(rand)= {f_rand_mean:.4f} +- {f_rand_std:.4f}  ({n_seeds} seeds)")
    print(f"  dF     = {advantage:+.4f}")
    if full_tomo and rho_full is not None:
        print(f"  F(full)= {f_full:.4f}  (gold standard)")
    print()
    print(f"  Z-only:  Delta_DTC={zonly_params['delta_dtc']:.4f}")
    print(f"  K* adds: XY_frac={kstar_direct_coh['xy_signal_fraction']:.4f}, "
          f"l1={coherence_kstar['l1_coherence']:.4f}, "
          f"neg={coherence_kstar['partial_transpose_negativity']:.4f}")
    if kstar_direct_coh['any_xy_nonzero']:
        print("  >> K* DETECTS XY COHERENCES BEYOND Z-ONLY PROTOCOL <<")
    print("=" * 70)

    return results


# ---------------------------------------------------------------------------
# Noiseless simulation mode
# ---------------------------------------------------------------------------

def run_simulation_only(n_qubits=4, phi=np.pi / 2, eps=0.1, J=1.0,
                        n_cycles=6, init_state='neel', n_shots=1000, seed=42,
                        n_seeds=5, noisy_sim=False,
                        depol_rate=0.02, readout_rate=0.02,
                        disorder_W=0.5, disorder_seed=0):
    """
    Run full analysis pipeline on simulated data.
    Validates the experiment design before using QPU time.

    Args:
        noisy_sim: if True, apply depolarizing + readout noise model
        depol_rate: per-qubit depolarizing probability (noisy mode)
        readout_rate: per-qubit readout error probability (noisy mode)
    """
    seeds = _make_seed_list(seed, n_seeds)
    bond_couplings = _generate_bond_couplings(n_qubits, J, disorder_W, disorder_seed)
    mode_str = f"noisy (depol={depol_rate}, readout={readout_rate})" if noisy_sim else "noiseless"

    print("=" * 70)
    print(f"  FLOQUET DTC -- SIMULATION ONLY ({mode_str})")
    print("=" * 70)
    print(f"  Params: phi={phi:.4f}, eps={eps}, J={J}, T={n_cycles}")
    print(f"  Disorder: W={disorder_W}, seed={disorder_seed}")
    bc_str = ', '.join(f"({i},{j}):{v:.3f}" for (i, j), v in sorted(bond_couplings.items()))
    print(f"  Bond J_ij: {bc_str}")
    print(f"  Init:   {init_state}, Shots: {n_shots}, Seeds: {seeds}")
    print()

    # Qubit ordering verification (C1)
    verify_qubit_ordering(verbose=True)
    print()

    # Ideal state
    rho_ideal = simulate_floquet_noiseless(n_qubits, phi, eps, J, n_cycles, init_state,
                                           bond_couplings=bond_couplings)

    # Compute all Pauli expectations from ideal state
    all_ops, all_labels, _ = all_pauli_operators(n_qubits)
    ideal_exps = {}
    for op, label in zip(all_ops, all_labels):
        ideal_exps[label] = np.trace(op @ rho_ideal).real

    # K* selection
    kstar_ops, kstar_labels, _ = select_kstar_paulis(n_qubits)

    # Noise function
    if noisy_sim:
        noisy_func = lambda exps, labels: simulate_noisy_expectations(
            exps, labels, n_shots, depol_rate, readout_rate)
        kstar_expectations = noisy_func(ideal_exps, kstar_labels)
    else:
        kstar_expectations = [shot_noise(ideal_exps[l], n_shots) for l in kstar_labels]
        noisy_func = None

    # Reconstruct K*
    rho_kstar = reconstruct_robust_mle(kstar_expectations, kstar_ops, n_qubits,
                                        n_shots=n_shots)
    f_kstar = state_fidelity(rho_ideal, rho_kstar)

    # Multi-seed random reconstruction
    f_rands = []
    for s in seeds:
        r_ops, r_labels, _ = select_random_paulis(n_qubits, len(kstar_ops), seed=s)
        if noisy_sim:
            re = noisy_func(ideal_exps, r_labels)
        else:
            re = [shot_noise(ideal_exps[l], n_shots) for l in r_labels]
        rho_r = reconstruct_robust_mle(re, r_ops, n_qubits, n_shots=n_shots)
        f_rands.append(state_fidelity(rho_ideal, rho_r))

    f_rand_mean = float(np.mean(f_rands))
    f_rand_std = float(np.std(f_rands))

    print(f"  F(K*)  = {f_kstar:.4f}")
    print(f"  F(rand)= {f_rand_mean:.4f} +- {f_rand_std:.4f}  ({n_seeds} seeds)")
    print(f"  dF     = {f_kstar - f_rand_mean:+.4f}")
    print()

    # Z-only order params
    z_exp, zz_exp = extract_zonly_from_full(kstar_expectations, kstar_labels, n_qubits)
    init_signs = _init_signs_for_state(init_state, n_qubits)
    zonly = compute_zonly_order_params(z_exp, zz_exp, n_qubits, init_signs, n_cycles)

    print(f"  Delta_DTC = {zonly['delta_dtc']:.4f}")
    print(f"  Delta_MBL = {zonly['delta_mbl']:.4f}")
    print()

    # K*-accessible coherences (model-independent, from expectations)
    kstar_direct = analyze_kstar_accessible_coherences(
        kstar_expectations, kstar_labels, n_qubits)
    print(f"  K*-direct XY signal fraction: {kstar_direct['xy_signal_fraction']:.4f}")
    print(f"  K*-direct any XY nonzero:     {kstar_direct['any_xy_nonzero']}")
    print()

    # Full coherence from reconstructed rho
    coh_kstar = analyze_coherences(rho_kstar, n_qubits)
    coh_ideal = analyze_coherences(rho_ideal, n_qubits)

    print(f"  l1 coherence (ideal): {coh_ideal['l1_coherence']:.4f}")
    print(f"  l1 coherence (K*):    {coh_kstar['l1_coherence']:.4f}")
    print(f"  PT negativity (ideal): {coh_ideal['partial_transpose_negativity']:.4f}")
    print(f"  PT negativity (K*):    {coh_kstar['partial_transpose_negativity']:.4f}")
    print()

    # K*-direct correlators
    print("  --- XY CORRELATORS (from K* expectations directly) ---")
    print(f"  XX: {kstar_direct['xx_correlators']}")
    print(f"  YY: {kstar_direct['yy_correlators']}")
    print(f"  XY: {kstar_direct['xy_correlators']}")
    print()

    # Weight-class breakdown
    print("  --- WEIGHT-CLASS EXPECTATION BREAKDOWN (ideal) ---")
    for w, ptypes in sorted(coh_ideal['weight_class_summary'].items()):
        for ptype, stats in ptypes.items():
            print(f"    w={w}, {ptype:8s}: mean|<P>|={stats['mean_abs']:.4f}, "
                  f"max={stats['max_abs']:.4f}, count={stats['count']}")

    # Stroboscopic trace
    strobe = compute_stroboscopic_trace(
        n_qubits, phi, eps, J, n_cycles, init_state, bond_couplings)
    print()
    print("  --- STROBOSCOPIC TRACE Delta(t) ---")
    for s in strobe:
        marker = " *" if s['delta_dtc'] > 0.5 else ""
        print(f"    t={s['t']:2d}: Delta={s['delta']:+.4f}, "
              f"Delta_DTC={s['delta_dtc']:+.4f}{marker}")

    # Hamming distribution
    hamming = compute_hamming_distribution(rho_ideal, init_state, n_qubits)
    print()
    print(f"  --- HAMMING DISTANCE (ideal, mean={hamming['mean_distance']:.2f}) ---")
    for h in hamming['distribution']:
        bar = '#' * int(h['probability'] * 40)
        print(f"    d={h['distance']}: {h['probability']:.4f} {bar}")

    print()
    print("  Simulation complete. Ready for hardware run.")
    print("=" * 70)


# ---------------------------------------------------------------------------
# Job recovery
# ---------------------------------------------------------------------------

def recover_from_jobs(api_token, recovery_file, instance=None):
    """
    Recover and process results from previously submitted jobs.
    Supports both single-point and sweep recovery files (auto-detected).
    """
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

    with open(recovery_file, 'r') as f:
        recovery_data = json.load(f)

    job_ids = recovery_data['job_ids']
    params = recovery_data['params']
    kstar_labels = recovery_data['kstar_labels']
    backend_name = recovery_data['backend']
    is_sweep = recovery_data.get('sweep', False)

    n_qubits = params['n_qubits']
    n_shots = params['n_shots']
    phi = params['phi']
    eps = params['eps']
    J = params['J']
    init_state = params['init_state']

    print(f"  Recovering {len(job_ids)} jobs from {backend_name}...")
    if is_sweep:
        print(f"  Mode: sweep (cycles={recovery_data['cycle_values']})")

    svc_kwargs = dict(channel="ibm_quantum_platform", token=api_token)
    if instance:
        svc_kwargs['instance'] = instance
    service = QiskitRuntimeService(**svc_kwargs)

    # Build circuit_labels in the same order as submission
    kstar_ops, _, _ = select_kstar_paulis(n_qubits)

    if is_sweep:
        cycle_values = recovery_data['cycle_values']
        n_seeds = recovery_data.get('n_seeds', 1)
        random_seeds = recovery_data.get('random_seeds', [42])
        rand_arms_meta = recovery_data.get('rand_arms', [])

        # Rebuild random operator sets from seeds
        rand_arms = []
        for arm_info in rand_arms_meta:
            s = arm_info['seed']
            r_ops, r_labels, _ = select_random_paulis(n_qubits, len(kstar_ops), seed=s)
            rand_arms.append((s, r_ops, r_labels))

        # Build ordered label list for expectation extraction
        circuit_labels = []
        for T in cycle_values:
            for label in kstar_labels:
                circuit_labels.append((T, 'kstar', label))
            for s, _, r_labels in rand_arms:
                for label in r_labels:
                    circuit_labels.append((T, f'random_s{s}', label))
    else:
        n_cycles = params['n_cycles']
        rand_labels = recovery_data['rand_labels']
        circuit_labels = [(0, 'kstar', l) for l in kstar_labels] + \
                         [(0, 'random', l) for l in rand_labels]

    # Retrieve all expectations
    all_expectations = []
    for job_id in job_ids:
        print(f"    Retrieving job {job_id}...")
        job = service.job(job_id)
        result = job.result()

        for j in range(len(result)):
            counts = None
            try:
                counts = result[j].data.meas.get_counts()
            except AttributeError:
                data_bin = result[j].data
                for attr_name in dir(data_bin):
                    if attr_name.startswith('_'):
                        continue
                    attr = getattr(data_bin, attr_name)
                    if hasattr(attr, 'get_counts'):
                        counts = attr.get_counts()
                        break

            if counts is None:
                all_expectations.append(0.0)
                continue

            idx = len(all_expectations)
            _, _, label = circuit_labels[idx]
            exp = counts_to_expectation(counts, label, n_qubits)
            all_expectations.append(exp)

    if is_sweep:
        # Process per-T with multi-seed random
        n_per_arm = len(kstar_labels)
        n_per_T = n_per_arm + n_seeds * n_per_arm
        init_signs = _init_signs_for_state(init_state, n_qubits)
        sweep_results = []

        for t_idx, T in enumerate(cycle_values):
            offset = t_idx * n_per_T
            ke = all_expectations[offset:offset + n_per_arm]

            rho_k = reconstruct_robust_mle(ke, kstar_ops, n_qubits, n_shots=n_shots)
            rho_ideal = simulate_floquet_noiseless(n_qubits, phi, eps, J, T, init_state)
            fk = state_fidelity(rho_ideal, rho_k)

            f_rands = []
            for seed_idx, (s, r_ops_s, r_labels_s) in enumerate(rand_arms):
                r_offset = offset + n_per_arm + seed_idx * n_per_arm
                re = all_expectations[r_offset:r_offset + n_per_arm]
                rho_r = reconstruct_robust_mle(re, r_ops_s, n_qubits, n_shots=n_shots)
                f_rands.append(state_fidelity(rho_ideal, rho_r))

            fr_mean = float(np.mean(f_rands))
            fr_std = float(np.std(f_rands))

            z_e, zz_e = extract_zonly_from_full(ke, kstar_labels, n_qubits)
            zp = compute_zonly_order_params(z_e, zz_e, n_qubits, init_signs, T)
            coh_k = analyze_coherences(rho_k, n_qubits)
            kstar_direct = analyze_kstar_accessible_coherences(ke, kstar_labels, n_qubits)

            row = {
                'T': T, 'f_kstar': float(fk),
                'f_rand_mean': fr_mean, 'f_rand_std': fr_std,
                'advantage_mean': float(fk - fr_mean),
                'delta_dtc': zp['delta_dtc'],
                'xy_signal_fraction': kstar_direct['xy_signal_fraction'],
            }
            sweep_results.append(row)
            print(f"  T={T}: F(K*)={fk:.4f}, F(rand)={fr_mean:.4f}+-{fr_std:.4f}, "
                  f"dF={fk-fr_mean:+.4f}")

        outfile = f"floquet_sweep_recovered_{time.strftime('%Y%m%d_%H%M%S')}.json"
        with open(outfile, 'w') as f:
            json.dump({'sweep': sweep_results, 'params': params,
                       'backend': backend_name}, f, indent=2)
        print(f"  Recovered sweep saved: {outfile}")
        return sweep_results

    else:
        # Single-point recovery (multi-seed if available)
        n_kstar = len(kstar_labels)
        kstar_expectations = all_expectations[:n_kstar]

        n_seeds = recovery_data.get('n_seeds', 1)
        random_seeds_list = recovery_data.get('random_seeds', [params.get('seed', 42)])
        rand_arms_meta = recovery_data.get('rand_arms', [])

        rho_kstar = reconstruct_robust_mle(
            kstar_expectations, kstar_ops, n_qubits, n_shots=n_shots)
        rho_ideal = simulate_floquet_noiseless(
            n_qubits, phi, eps, J, n_cycles, init_state)
        f_kstar = state_fidelity(rho_ideal, rho_kstar)

        if rand_arms_meta:
            # Multi-seed random
            f_rands = []
            pos = n_kstar
            for arm_info in rand_arms_meta:
                s = arm_info['seed']
                r_ops, r_labels, _ = select_random_paulis(
                    n_qubits, len(kstar_ops), seed=s)
                re = all_expectations[pos:pos + len(r_labels)]
                pos += len(r_labels)
                rho_r = reconstruct_robust_mle(re, r_ops, n_qubits, n_shots=n_shots)
                f_rands.append(state_fidelity(rho_ideal, rho_r))

            f_rand_mean = float(np.mean(f_rands))
            f_rand_std = float(np.std(f_rands))
            print(f"  F(K*) = {f_kstar:.4f}, F(rand) = {f_rand_mean:.4f}+-{f_rand_std:.4f}, "
                  f"dF = {f_kstar - f_rand_mean:+.4f}")
            return {'f_kstar': f_kstar, 'f_rand_mean': f_rand_mean,
                    'f_rand_std': f_rand_std, 'advantage': f_kstar - f_rand_mean}
        else:
            # Legacy single-seed recovery
            rand_labels = recovery_data['rand_labels']
            rand_expectations = all_expectations[n_kstar:]
            rand_ops, _, _ = select_random_paulis(n_qubits, len(kstar_ops), seed=42)
            rho_rand = reconstruct_robust_mle(
                rand_expectations, rand_ops, n_qubits, n_shots=n_shots)
            f_rand = state_fidelity(rho_ideal, rho_rand)
            print(f"  F(K*) = {f_kstar:.4f}, F(rand) = {f_rand:.4f}, "
                  f"dF = {f_kstar - f_rand:+.4f}")
            return {'f_kstar': f_kstar, 'f_rand': f_rand,
                    'advantage': f_kstar - f_rand}


# ---------------------------------------------------------------------------
# Multi-cycle sweep (main figure data)
# ---------------------------------------------------------------------------

def run_multicycle_sweep(api_token=None, n_qubits=4, n_shots=1000,
                         backend_name=None, initial_qubits=None,
                         phi=np.pi / 2, eps=0.1, J=1.0,
                         init_state='neel', seed=42, n_seeds=5,
                         sim_only=False, noisy_sim=False,
                         depol_rate=0.02, readout_rate=0.02,
                         cycle_values=None, grouped=False,
                         disorder_W=0.5, disorder_seed=0,
                         instance=None):
    """
    Sweep over Floquet cycle counts T=2,4,6,8, collecting K* vs random
    fidelity, Z-only order parameters, and coherence metrics at each T.

    Produces the main figure data: dF(T), Delta_DTC(T), l1_coherence(T).
    Uses multi-seed random baseline for error bars.
    """
    if cycle_values is None:
        cycle_values = [2, 4, 6, 8]

    seeds = _make_seed_list(seed, n_seeds)
    bond_couplings = _generate_bond_couplings(n_qubits, J, disorder_W, disorder_seed)
    mode_str = 'simulation'
    if sim_only and noisy_sim:
        mode_str = f'noisy sim (depol={depol_rate}, readout={readout_rate})'
    elif not sim_only:
        mode_str = 'hardware'

    print("=" * 70)
    print("  FLOQUET DTC -- MULTI-CYCLE SWEEP")
    print("=" * 70)
    print(f"  Cycles:  {cycle_values}")
    print(f"  Params:  phi={phi:.4f}, eps={eps}, J={J}")
    print(f"  Disorder: W={disorder_W}, seed={disorder_seed}")
    bc_str = ', '.join(f"({i},{j}):{v:.3f}" for (i, j), v in sorted(bond_couplings.items()))
    print(f"  Bond J_ij: {bc_str}")
    print(f"  Seeds:   {seeds}")
    print(f"  Mode:    {mode_str}")
    print()

    kstar_ops, kstar_labels, _ = select_kstar_paulis(n_qubits)
    all_ops, all_labels, _ = all_pauli_operators(n_qubits)
    init_signs = _init_signs_for_state(init_state, n_qubits)

    sweep_results = []

    if sim_only:
        if noisy_sim:
            noisy_func = lambda exps, labels: simulate_noisy_expectations(
                exps, labels, n_shots, depol_rate, readout_rate)
        else:
            noisy_func = None

        for T in cycle_values:
            print(f"  --- T={T} ---")
            rho_ideal = simulate_floquet_noiseless(n_qubits, phi, eps, J, T, init_state,
                                                   bond_couplings=bond_couplings)

            ideal_exps = {l: np.trace(o @ rho_ideal).real
                          for o, l in zip(all_ops, all_labels)}

            if noisy_func:
                ke = noisy_func(ideal_exps, kstar_labels)
            else:
                ke = [shot_noise(ideal_exps[l], n_shots) for l in kstar_labels]

            rho_k = reconstruct_robust_mle(ke, kstar_ops, n_qubits, n_shots=n_shots)
            fk = state_fidelity(rho_ideal, rho_k)

            # Multi-seed random
            f_rands = []
            for s in seeds:
                r_ops, r_labels, _ = select_random_paulis(n_qubits, len(kstar_ops), seed=s)
                if noisy_func:
                    re = noisy_func(ideal_exps, r_labels)
                else:
                    re = [shot_noise(ideal_exps[l], n_shots) for l in r_labels]
                rho_r = reconstruct_robust_mle(re, r_ops, n_qubits, n_shots=n_shots)
                f_rands.append(state_fidelity(rho_ideal, rho_r))

            fr_mean = float(np.mean(f_rands))
            fr_std = float(np.std(f_rands))

            z_e, zz_e = extract_zonly_from_full(ke, kstar_labels, n_qubits)
            zp = compute_zonly_order_params(z_e, zz_e, n_qubits, init_signs, T)

            coh_k = analyze_coherences(rho_k, n_qubits)
            coh_ideal = analyze_coherences(rho_ideal, n_qubits)
            kstar_direct = analyze_kstar_accessible_coherences(ke, kstar_labels, n_qubits)

            # Stroboscopic trace and Hamming for this T
            strobe = compute_stroboscopic_trace(
                n_qubits, phi, eps, J, T, init_state, bond_couplings)
            hamming = compute_hamming_distribution(rho_ideal, init_state, n_qubits)
            hamming_k = compute_hamming_distribution(rho_k, init_state, n_qubits)

            row = {
                'T': T,
                'f_kstar': float(fk),
                'f_rand_mean': fr_mean,
                'f_rand_std': fr_std,
                'f_rand_per_seed': {str(s): float(f) for s, f in zip(seeds, f_rands)},
                'advantage_mean': float(fk - fr_mean),
                'delta_dtc': zp['delta_dtc'],
                'delta_mbl': zp['delta_mbl'],
                'l1_coherence_ideal': coh_ideal['l1_coherence'],
                'l1_coherence_kstar': coh_k['l1_coherence'],
                'negativity_ideal': coh_ideal['partial_transpose_negativity'],
                'negativity_kstar': coh_k['partial_transpose_negativity'],
                'xy_signal_fraction': kstar_direct['xy_signal_fraction'],
                'xx_corr_mean': float(np.mean([abs(v) for v in coh_k['xx_correlators'].values()])),
                'yy_corr_mean': float(np.mean([abs(v) for v in coh_k['yy_correlators'].values()])),
                'stroboscopic_trace': strobe,
                'hamming_ideal': hamming,
                'hamming_kstar': hamming_k,
            }
            sweep_results.append(row)

            print(f"    F(K*)={fk:.4f}, F(rand)={fr_mean:.4f}+-{fr_std:.4f}, "
                  f"dF={fk-fr_mean:+.4f}, DTC={zp['delta_dtc']:+.3f}")

    else:
        # Hardware mode: build all circuits, submit in one session
        if not api_token:
            print("ERROR: API token required for hardware sweep")
            return

        from qiskit import QuantumCircuit
        from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
        from qiskit.transpiler import preset_passmanagers

        print("  Connecting to IBM Quantum...")
        svc_kwargs = dict(channel="ibm_quantum_platform", token=api_token)
        if instance:
            svc_kwargs['instance'] = instance
        service = QiskitRuntimeService(**svc_kwargs)
        if backend_name:
            backend = service.backend(backend_name)
        else:
            backends = service.backends(
                filters=lambda b: b.num_qubits >= n_qubits and b.status().operational)
            backend = min(backends, key=lambda b: b.status().pending_jobs)
        print(f"  Backend: {backend.name}")

        # Build random arm label sets for each seed
        rand_arms = []
        for s in seeds:
            r_ops, r_labels, _ = select_random_paulis(n_qubits, len(kstar_ops), seed=s)
            rand_arms.append((s, r_ops, r_labels))

        # Build all circuits for all T values
        all_circuits = []
        circuit_metadata = []  # (T, arm, label_or_basis)
        # For grouped mode, we also need per-T metadata for extraction
        sweep_grouped_meta = []  # per-circuit: (T, arm, basis, [covered_labels])

        if grouped:
            for T in cycle_values:
                label_sets = [('kstar', list(kstar_labels))]
                for s, _, r_labels in rand_arms:
                    label_sets.append((f'random_s{s}', list(r_labels)))
                for arm_name, labels in label_sets:
                    cover = greedy_basis_cover(labels, n_qubits)
                    for basis, covered in cover:
                        qc = floquet_circuit_with_measurement(
                            n_qubits, phi, eps, J, T, basis, init_state,
                            bond_couplings=bond_couplings)
                        qc.name = f"T{T}_{arm_name}_{basis}"
                        all_circuits.append(qc)
                        circuit_metadata.append((T, arm_name, basis))
                        sweep_grouped_meta.append((T, arm_name, basis, covered))
        else:
            for T in cycle_values:
                for label in kstar_labels:
                    qc = floquet_circuit_with_measurement(
                        n_qubits, phi, eps, J, T, label, init_state,
                        bond_couplings=bond_couplings)
                    qc.name = f"T{T}_kstar_{label}"
                    all_circuits.append(qc)
                    circuit_metadata.append((T, 'kstar', label))

                for s, r_ops, r_labels in rand_arms:
                    for label in r_labels:
                        qc = floquet_circuit_with_measurement(
                            n_qubits, phi, eps, J, T, label, init_state,
                            bond_couplings=bond_couplings)
                        qc.name = f"T{T}_rand_s{s}_{label}"
                        all_circuits.append(qc)
                        circuit_metadata.append((T, f'random_s{s}', label))

        n_per_arm = len(kstar_labels)
        total = len(all_circuits)
        if grouped:
            print(f"  GROUPED: {total} circuits (vs {len(cycle_values) * (n_per_arm + n_seeds * n_per_arm)} ungrouped)")
        else:
            n_per_T = n_per_arm + n_seeds * n_per_arm
            print(f"  Total circuits: {total} ({len(cycle_values)} x (K*={n_per_arm} + rand={n_seeds}x{n_per_arm}))")

        # Transpile
        print("  Transpiling...")
        pm_kwargs = {'optimization_level': 1, 'backend': backend}
        if initial_qubits:
            pm_kwargs['initial_layout'] = initial_qubits
        pm = preset_passmanagers.generate_preset_pass_manager(**pm_kwargs)
        transpiled = pm.run(all_circuits)

        # Submit
        print("  Submitting...")
        batch_size = 300 if grouped else 100
        all_hw_results = []
        job_ids = []
        sampler = SamplerV2(mode=backend)

        for i in range(0, len(transpiled), batch_size):
            batch = transpiled[i:i + batch_size]
            batch_pubs = [(circ,) for circ in batch]
            job = sampler.run(batch_pubs)
            job_ids.append(job.job_id())
            print(f"    Batch {i // batch_size + 1}: {len(batch)} circuits, job={job.job_id()}")
            t0 = time.time()
            result = job.result()
            print(f"    -> {time.time() - t0:.1f}s")
            all_hw_results.append((i, result))

        # Save recovery
        recovery = {
            'job_ids': job_ids,
            'backend': backend.name,
            'sweep': True,
            'cycle_values': cycle_values,
            'n_seeds': n_seeds,
            'random_seeds': seeds,
            'params': {'phi': phi, 'eps': eps, 'J': J, 'init_state': init_state,
                       'n_qubits': n_qubits, 'n_shots': n_shots},
            'kstar_labels': list(kstar_labels),
            'rand_arms': [{'seed': s, 'labels': list(rl)} for s, _, rl in rand_arms],
        }
        rf = f"floquet_sweep_jobs_{backend.name}_{time.strftime('%Y%m%d_%H%M%S')}.json"
        with open(rf, 'w') as f:
            json.dump(recovery, f, indent=2)

        # Extract expectations
        all_counts = []
        for batch_start, result in all_hw_results:
            for j in range(len(result)):
                pub_result = result[j]
                counts = None
                try:
                    counts = pub_result.data.meas.get_counts()
                except AttributeError:
                    data_bin = pub_result.data
                    for attr_name in dir(data_bin):
                        if attr_name.startswith('_'):
                            continue
                        attr = getattr(data_bin, attr_name)
                        if hasattr(attr, 'get_counts'):
                            counts = attr.get_counts()
                            break
                all_counts.append(counts if counts else {})

        if grouped:
            # Grouped extraction: fan out per-operator expectations per T
            # Build per-T arm expectations dicts
            per_T_exps = {T: defaultdict(dict) for T in cycle_values}
            for ci, (T, arm, basis, covered) in enumerate(sweep_grouped_meta):
                for label in covered:
                    per_T_exps[T][arm][label] = counts_to_expectation(
                        all_counts[ci], label, n_qubits)
        else:
            # Ungrouped: convert counts to expectations using metadata labels
            all_expectations = []
            for ci, counts in enumerate(all_counts):
                if not counts:
                    all_expectations.append(0.0)
                    continue
                _, _, label = circuit_metadata[ci]
                all_expectations.append(counts_to_expectation(counts, label, n_qubits))

        # Process per-T (multi-seed random)
        n_per_arm = len(kstar_labels)

        for t_idx, T in enumerate(cycle_values):
            if grouped:
                t_exps = per_T_exps[T]
                ke = [t_exps.get('kstar', {}).get(l, 0.0) for l in kstar_labels]
            else:
                n_per_T = n_per_arm + n_seeds * n_per_arm
                offset = t_idx * n_per_T
                ke = all_expectations[offset:offset + n_per_arm]

            rho_k = reconstruct_robust_mle(ke, kstar_ops, n_qubits, n_shots=n_shots)
            rho_ideal = simulate_floquet_noiseless(n_qubits, phi, eps, J, T, init_state,
                                                   bond_couplings=bond_couplings)
            fk = state_fidelity(rho_ideal, rho_k)

            # Multi-seed random fidelities
            f_rands = []
            for seed_idx, (s, r_ops_s, r_labels_s) in enumerate(rand_arms):
                if grouped:
                    arm_key = f'random_s{s}'
                    re = [t_exps.get(arm_key, {}).get(l, 0.0) for l in r_labels_s]
                else:
                    r_offset = offset + n_per_arm + seed_idx * n_per_arm
                    re = all_expectations[r_offset:r_offset + n_per_arm]
                rho_r = reconstruct_robust_mle(re, r_ops_s, n_qubits, n_shots=n_shots)
                f_rands.append(state_fidelity(rho_ideal, rho_r))

            fr_mean = float(np.mean(f_rands))
            fr_std = float(np.std(f_rands))

            z_e, zz_e = extract_zonly_from_full(ke, kstar_labels, n_qubits)
            zp = compute_zonly_order_params(z_e, zz_e, n_qubits, init_signs, T)

            coh_k = analyze_coherences(rho_k, n_qubits)
            coh_ideal = analyze_coherences(rho_ideal, n_qubits)
            kstar_direct = analyze_kstar_accessible_coherences(ke, kstar_labels, n_qubits)

            strobe = compute_stroboscopic_trace(
                n_qubits, phi, eps, J, T, init_state, bond_couplings)
            hamming_hw = compute_hamming_distribution(rho_k, init_state, n_qubits)
            hamming_id = compute_hamming_distribution(rho_ideal, init_state, n_qubits)

            row = {
                'T': T,
                'f_kstar': float(fk),
                'f_rand_mean': fr_mean,
                'f_rand_std': fr_std,
                'f_rand_per_seed': {str(s): float(f) for s, f in zip(seeds, f_rands)},
                'advantage_mean': float(fk - fr_mean),
                'delta_dtc': zp['delta_dtc'],
                'delta_mbl': zp['delta_mbl'],
                'l1_coherence_ideal': coh_ideal['l1_coherence'],
                'l1_coherence_kstar': coh_k['l1_coherence'],
                'negativity_ideal': coh_ideal['partial_transpose_negativity'],
                'negativity_kstar': coh_k['partial_transpose_negativity'],
                'xy_signal_fraction': kstar_direct['xy_signal_fraction'],
                'xx_corr_mean': float(np.mean([abs(v) for v in coh_k['xx_correlators'].values()])),
                'yy_corr_mean': float(np.mean([abs(v) for v in coh_k['yy_correlators'].values()])),
                'kstar_expectations': [float(e) for e in ke],
                'stroboscopic_trace': strobe,
                'hamming_kstar': hamming_hw,
                'hamming_ideal': hamming_id,
            }
            sweep_results.append(row)

            print(f"  T={T}: F(K*)={fk:.4f}, F(rand)={fr_mean:.4f}+-{fr_std:.4f}, "
                  f"dF={fk-fr_mean:+.4f}, DTC={zp['delta_dtc']:+.3f}")

    # Save sweep results
    output = {
        'experiment': 'floquet_dtc_multicycle_sweep',
        'mode': 'simulation' if sim_only else 'hardware',
        'params': {'phi': float(phi), 'eps': float(eps), 'J': float(J),
                   'init_state': init_state, 'n_shots': n_shots,
                   'disorder_W': float(disorder_W), 'disorder_seed': int(disorder_seed),
                   'bond_couplings': {f"{k[0]},{k[1]}": float(v)
                                      for k, v in bond_couplings.items()}},
        'sweep': sweep_results,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    }
    if not sim_only:
        output['backend'] = backend.name
        output['job_ids'] = job_ids

    outfile = f"floquet_sweep_results_{time.strftime('%Y%m%d_%H%M%S')}.json"
    with open(outfile, 'w') as f:
        json.dump(output, f, indent=2)

    # Print summary table
    print()
    print("  === SWEEP SUMMARY ===")
    print(f"  {'T':>3s}  {'F(K*)':>7s}  {'F(rand)':>12s}  {'dF':>7s}  {'DTC':>6s}  {'l1_coh':>7s}  {'XY_frac':>7s}")
    print("  " + "-" * 65)
    for r in sweep_results:
        fr_str = f"{r['f_rand_mean']:.4f}+-{r['f_rand_std']:.3f}"
        print(f"  {r['T']:3d}  {r['f_kstar']:7.4f}  {fr_str:>12s}  {r['advantage_mean']:+7.4f}  "
              f"{r['delta_dtc']:+6.3f}  {r['l1_coherence_kstar']:7.3f}  {r['xy_signal_fraction']:7.4f}")

    print()
    print(f"  Results saved: {outfile}")
    print("=" * 70)

    return sweep_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Floquet DTC K* Certification on IBM Hardware'
    )
    parser.add_argument('api_token', nargs='?', default=None,
                        help='IBM Quantum API token')
    parser.add_argument('--backend', type=str, default=None,
                        help='IBM backend name')
    parser.add_argument('--instance', type=str, default=None,
                        help='IBM Cloud instance CRN (required for new platform accounts)')
    parser.add_argument('--qubits', type=str, default=None,
                        help='Physical qubit indices (comma-separated)')
    parser.add_argument('--phi', type=float, default=np.pi / 2,
                        help='Floquet kick angle (default: pi/2)')
    parser.add_argument('--eps', type=float, default=0.1,
                        help='Heisenberg anisotropy (default: 0.1)')
    parser.add_argument('--J', type=float, default=1.0,
                        help='Base coupling strength (default: 1.0, matches BASQ)')
    parser.add_argument('--cycles', type=int, default=6,
                        help='Number of Floquet cycles (default: 6)')
    parser.add_argument('--n-shots', type=int, default=1000,
                        help='Shots per circuit (default: 1000)')
    parser.add_argument('--init-state', type=str, default='neel',
                        choices=['neel', 'afm', 'polarized'],
                        help='Initial state (default: neel)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for random Paulis')
    parser.add_argument('--n-seeds', type=int, default=5,
                        help='Number of random seeds for multi-seed baseline (default: 5)')
    parser.add_argument('--recover', type=str, default=None,
                        help='Recovery JSON file from previous run')
    parser.add_argument('--sim-only', action='store_true',
                        help='Noiseless simulation only (no IBM)')
    parser.add_argument('--sweep', action='store_true',
                        help='Multi-cycle sweep (main figure data)')
    parser.add_argument('--cycle-values', type=str, default=None,
                        help='Comma-separated cycle counts for sweep (default: 2,4,6,8)')
    parser.add_argument('--noisy-sim', action='store_true',
                        help='Add depolarizing + readout noise to simulation')
    parser.add_argument('--full-tomo', action='store_true',
                        help='Include full-tomography gold standard arm')
    parser.add_argument('--grouped', action='store_true',
                        help='Use basis grouping to reduce circuit count (29 vs 137 for K*)')
    parser.add_argument('--depol-rate', type=float, default=0.02,
                        help='Depolarizing error rate for noisy sim (default: 0.02)')
    parser.add_argument('--readout-rate', type=float, default=0.02,
                        help='Readout error rate for noisy sim (default: 0.02)')
    parser.add_argument('--disorder', type=float, default=0.5,
                        help='Disorder width W: J_ij = J*(1+delta), delta~U[-W/2,W/2] '
                             '(default: 0.5, matches BASQ; 0 for uniform)')
    parser.add_argument('--disorder-seed', type=int, default=0,
                        help='RNG seed for disorder realization (default: 0)')

    args = parser.parse_args()

    initial_qubits = ([int(q) for q in args.qubits.split(',')]
                      if args.qubits else None)

    cycle_values = None
    if args.cycle_values:
        cycle_values = [int(x) for x in args.cycle_values.split(',')]

    if args.sweep:
        run_multicycle_sweep(
            api_token=args.api_token,
            n_qubits=4, n_shots=args.n_shots,
            backend_name=args.backend,
            initial_qubits=initial_qubits,
            phi=args.phi, eps=args.eps, J=args.J,
            init_state=args.init_state, seed=args.seed,
            n_seeds=args.n_seeds,
            sim_only=args.sim_only,
            noisy_sim=args.noisy_sim,
            depol_rate=args.depol_rate,
            readout_rate=args.readout_rate,
            cycle_values=cycle_values,
            grouped=args.grouped,
            disorder_W=args.disorder,
            disorder_seed=args.disorder_seed,
            instance=args.instance,
        )
        return

    if args.sim_only:
        run_simulation_only(
            n_qubits=4, phi=args.phi, eps=args.eps, J=args.J,
            n_cycles=args.cycles, init_state=args.init_state,
            n_shots=args.n_shots, seed=args.seed,
            n_seeds=args.n_seeds,
            noisy_sim=args.noisy_sim,
            depol_rate=args.depol_rate,
            readout_rate=args.readout_rate,
            disorder_W=args.disorder,
            disorder_seed=args.disorder_seed,
        )
        return

    if args.recover:
        if not args.api_token:
            print("ERROR: API token required for recovery")
            sys.exit(1)
        recover_from_jobs(args.api_token, args.recover, instance=args.instance)
        return

    if not args.api_token:
        print("ERROR: API token required (or use --sim-only)")
        sys.exit(1)

    run_floquet_dtc_test(
        api_token=args.api_token,
        n_qubits=4,
        n_shots=args.n_shots,
        backend_name=args.backend,
        initial_qubits=initial_qubits,
        phi=args.phi,
        eps=args.eps,
        J=args.J,
        n_cycles=args.cycles,
        init_state=args.init_state,
        seed=args.seed,
        n_seeds=args.n_seeds,
        full_tomo=args.full_tomo,
        grouped=args.grouped,
        disorder_W=args.disorder,
        disorder_seed=args.disorder_seed,
        instance=args.instance,
    )


if __name__ == '__main__':
    main()
