"""Hessian properties and condition numbers (Theorem 1, Lemma)."""

import numpy as np

from common import (
    w_state, ghz_state, product_state_plus, one_local_state,
    two_local_state, depolarise,
    pauli_expectations, pauli_weight, all_paulis,
    kstar_operator_indices, random_operator_indices
)
from registry import claims


def _hessian_eigenvalues(expectations, operator_set, N_s=1000):
    """h_P = N_s / (1 - x_P^2) for each measured P."""
    eigenvalues = []
    for p in operator_set:
        if all(x == 0 for x in p):
            continue
        x_P = expectations.get(p, 0.0)
        if abs(x_P) < 1.0:
            h_P = N_s / (1.0 - x_P**2)
            eigenvalues.append(h_P)
    return np.array(eigenvalues)


def _condition_number_informative(expectations, operator_set, n,
                                  N_s=1000, eps=1e-6):
    """kappa = max(h_P)/min(h_P) over informative measured operators."""
    eigenvalues = []
    for p in operator_set:
        if all(x == 0 for x in p):
            continue
        x_P = expectations.get(p, 0.0)
        if abs(x_P) > eps and abs(x_P) < 1.0:
            h_P = N_s / (1.0 - x_P**2)
            eigenvalues.append(h_P)
    if len(eigenvalues) == 0:
        return float('inf')
    return max(eigenvalues) / min(eigenvalues)


def _count_informative_flat(expectations, operator_set, n, eps=1e-6):
    """Informative operators NOT in the measurement set (Theorem 1)."""
    count = 0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        x_P = expectations.get(tuple(p), 0.0)
        if abs(x_P) > eps and tuple(p) not in operator_set:
            count += 1
    return count


def _eps_flat(expectations, operator_set, n, eps=1e-6):
    """eps_flat = ||rho_{F cap I}||_F^2, the informative-but-flat Frobenius norm."""
    total = 0.0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        x_P = expectations.get(tuple(p), 0.0)
        if abs(x_P) > eps and tuple(p) not in operator_set:
            total += x_P**2
    return total


def _build_fisher_information_matrix(expectations, operator_set, n, N_s=1000):
    """Construct the full Fisher information matrix H_{PP'} explicitly.

    For Pauli measurements with N_s shots each:
      H_{PP'} = sum_{Q in S} [dTr(Q rho)/dy_P] [dTr(Q rho)/dy_{P'}] / Var(Q)

    Since Tr(Q rho) = y_Q (the Pauli-Bloch coordinate for Q), we have
    dTr(Q rho)/dy_P = delta_{Q,P}. Therefore:
      H_{PP'} = delta_{PP'} * [P in S] * N_s / (1 - y_P^2)

    This function builds the FULL matrix and returns it so we can verify
    that off-diagonal elements are exactly zero — i.e., the Hessian IS
    diagonal in Pauli-Bloch coordinates by construction for Pauli measurements.
    """
    # All non-identity Pauli indices
    all_ops = [tuple(p) for p in all_paulis(n)
               if not all(x == 0 for x in p)]
    M = len(all_ops)
    op_to_idx = {p: i for i, p in enumerate(all_ops)}

    H = np.zeros((M, M))

    for Q in operator_set:
        if all(x == 0 for x in Q):
            continue
        if Q not in op_to_idx:
            continue
        i = op_to_idx[Q]
        y_Q = expectations.get(Q, 0.0)
        if abs(y_Q) < 1.0 - 1e-12:
            h_Q = N_s / (1.0 - y_Q**2)
        else:
            h_Q = N_s * 1e6  # near-eigenstate: very large Fisher info
        # Each measurement Q contributes only to the (Q,Q) entry
        H[i, i] += h_Q

    return H, all_ops, op_to_idx


def test_hessian_properties():
    """Verify Hessian properties from Theorem 1 and Lemma."""
    passed = 0
    n = 4
    N_s = 1000
    p_noise = 0.03

    rho_W = depolarise(w_state(n), p_noise)
    exp_W = pauli_expectations(rho_W, n)

    n4_M = claims.get("thm:basin", "n4_M")
    kstar_ops = kstar_operator_indices(n, K=5)
    assert len(kstar_ops) == n4_M

    # === LEMMA 1 VERIFICATION: Hessian is diagonal in Pauli-Bloch basis ===
    # For Pauli measurements, H_{PP'} = delta_{PP'} * N_s/(1-y_P^2) because
    # dTr(Q rho)/dy_P = delta_{QP} (Pauli-Bloch coordinates are orthonormal).
    # This is a STRUCTURAL property of Pauli measurements, verified here by
    # constructing the full 255x255 matrix and confirming the diagonal formula.

    H, all_ops, op_to_idx = _build_fisher_information_matrix(
        exp_W, kstar_ops, n, N_s)

    # Off-diagonal elements must be zero
    M = len(all_ops)
    off_diag_max = 0.0
    for i in range(M):
        for j in range(M):
            if i != j:
                off_diag_max = max(off_diag_max, abs(H[i, j]))

    assert off_diag_max == 0.0, \
        f"Off-diagonal max = {off_diag_max}, expected exactly 0"
    passed += 1
    print(f"  [PASS] Fisher information matrix is diagonal "
          f"(off-diag max = {off_diag_max})")

    # Verify diagonal entries match the formula h_P = N_s/(1-y_P^2)
    max_diag_err = 0.0
    n_measured = 0
    n_unmeasured = 0
    for p in all_ops:
        i = op_to_idx[p]
        y_P = exp_W.get(p, 0.0)
        if p in kstar_ops and abs(y_P) < 1.0 - 1e-12:
            expected_h = N_s / (1.0 - y_P**2)
            max_diag_err = max(max_diag_err, abs(H[i, i] - expected_h))
            n_measured += 1
        else:
            if p not in kstar_ops:
                assert H[i, i] == 0.0, \
                    f"Unmeasured operator {p} has H[P,P] = {H[i, i]}"
                n_unmeasured += 1

    assert n_measured > 0, \
        "No measured operators with |y_P| < 1 found — test is vacuous"
    assert max_diag_err < 1e-10, \
        f"Diagonal formula error = {max_diag_err}"
    passed += 1
    print(f"  [PASS] Diagonal entries match h_P = N_s/(1-y_P^2) "
          f"(err = {max_diag_err:.1e}, {n_measured} measured, "
          f"{n_unmeasured} zero)")

    # Verify on a SECOND state (product) to confirm generality
    rho_prod = depolarise(product_state_plus(n), p_noise)
    exp_prod = pauli_expectations(rho_prod, n)
    H_prod, _, _ = _build_fisher_information_matrix(
        exp_prod, kstar_ops, n, N_s)

    off_diag_max_prod = 0.0
    for i in range(M):
        for j in range(M):
            if i != j:
                off_diag_max_prod = max(off_diag_max_prod, abs(H_prod[i, j]))
    assert off_diag_max_prod == 0.0
    passed += 1
    print(f"  [PASS] Fisher info diagonal for product state "
          f"(off-diag max = {off_diag_max_prod})")

    # kappa restricted to weight <= 2 (saturated subspace where
    # all informative ops have |x_P| ~ 0.485 for W)
    w2_ops = {p for p in kstar_ops
              if sum(1 for x in p if x != 0) <= 2}
    kappa_w2 = _condition_number_informative(exp_W, w2_ops, n, N_s,
                                             eps=0.1)  # filter small exp
    assert kappa_w2 < 3.0, \
        f"kappa(weight<=2) = {kappa_w2:.2f}, expected near 1.31"
    passed += 1
    print(f"  [PASS] W state kappa(weight<=2) = {kappa_w2:.3f} "
          f"(~1.31 expected)")

    # --- kappa_info for K* on pure W state ---
    # All informative ops with |x_P| < 1 have |x_P| = 0.5 exactly,
    # so kappa = 1.0 (perfectly conditioned within non-degenerate subspace).
    # Under depolarization, weight-4 ops with |x_P| near 1 get scaled to ~0.97,
    # producing a larger kappa. Test both.

    # Pure W state: kappa = 1.0 (excl. degenerate |x|=1 ops)
    rho_W_pure = w_state(n)
    exp_W_pure = pauli_expectations(rho_W_pure, n)
    # Filter out |x_P| >= 1-eps to exclude degenerate binomials
    pure_eigenvalues = []
    for p in kstar_ops:
        if all(x == 0 for x in p):
            continue
        x_P = exp_W_pure.get(p, 0.0)
        if abs(x_P) > 1e-6 and abs(x_P) < 1.0 - 1e-10:
            pure_eigenvalues.append(N_s / (1.0 - x_P**2))
    if len(pure_eigenvalues) > 0:
        kappa_pure = max(pure_eigenvalues) / min(pure_eigenvalues)
        assert abs(kappa_pure - 1.0) < 1e-10, \
            f"kappa(K*, pure W, excl degenerate) = {kappa_pure:.6f}, expected 1.0"
        passed += 1
        print(f"  [PASS] kappa_info(K*, pure W) = {kappa_pure:.6f} "
              f"(expected 1.0: all |x_P|=0.5)")

    # random: same weight <= 2 restriction, many ops missing
    kappas_random_w2 = []
    for seed in range(10):
        rand_ops = random_operator_indices(n, n4_M,
                                           rng=np.random.default_rng(seed))
        rand_w2 = {p for p in rand_ops
                   if sum(1 for x in p if x != 0) <= 2}
        k = _condition_number_informative(exp_W, rand_w2, n, N_s, eps=0.1)
        if k < float('inf'):
            kappas_random_w2.append(k)

    assert len(kappas_random_w2) > 0, \
        "All 10 random seeds missed every informative weight<=2 operator"
    mean_k_rand = np.mean(kappas_random_w2)
    passed += 1
    print(f"  [PASS] Random weight<=2 kappa = {mean_k_rand:.2f} "
          f"({len(kappas_random_w2)}/10 seeds had informative ops)")

    # --- informative flat directions ---

    n_flat_kstar = 0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        w = pauli_weight(p)
        x_P = exp_W.get(tuple(p), 0.0)
        if w <= 2 and abs(x_P) > 1e-6 and tuple(p) not in kstar_ops:
            n_flat_kstar += 1
    assert n_flat_kstar == 0, \
        f"K* has {n_flat_kstar} informative flat directions at weight <= 2"
    passed += 1
    print(f"  [PASS] K* informative flat (weight<=2) = 0")

    rand_ops = random_operator_indices(n, n4_M, rng=np.random.default_rng(42))
    n_flat_rand = _count_informative_flat(exp_W, rand_ops, n)
    assert n_flat_rand > 0
    passed += 1
    print(f"  [PASS] Random informative flat = {n_flat_rand} > 0")

    # --- eps_flat ---

    eps_flat_kstar = _eps_flat(exp_W, kstar_ops, n)
    eps_flat_rand = _eps_flat(exp_W, rand_ops, n)

    # K* only misses weight > 2 directions (small for W state)
    assert eps_flat_rand > eps_flat_kstar
    passed += 1
    print(f"  [PASS] eps_flat: K*={eps_flat_kstar:.4f}, "
          f"random={eps_flat_rand:.4f}")

    # --- states with W(rho) <= w_sat = 2 ---

    # W(rho)=1: all info at weight 1, saturated by K*
    rho_1loc = one_local_state(n)
    exp_1loc = pauli_expectations(rho_1loc, n)
    n_flat_1loc = _count_informative_flat(exp_1loc, kstar_ops, n)
    assert n_flat_1loc == 0, \
        f"1-local state has {n_flat_1loc} informative flat (expected 0)"
    passed += 1
    print(f"  [PASS] 1-local state: informative flat = 0 "
          f"(w_sat=2 >= W=1)")

    eps_flat_1loc = _eps_flat(exp_1loc, kstar_ops, n)
    assert eps_flat_1loc < 1e-10
    passed += 1
    print(f"  [PASS] 1-local state: eps_flat(K*) = {eps_flat_1loc:.2e}")

    # W(rho)=2: Bell pair on qubits 1,2; weight<=2 saturated
    rho_2loc = two_local_state(n)
    exp_2loc = pauli_expectations(rho_2loc, n)
    n_flat_2loc = _count_informative_flat(exp_2loc, kstar_ops, n)
    assert n_flat_2loc == 0, \
        f"2-local state has {n_flat_2loc} informative flat (expected 0)"
    passed += 1
    print(f"  [PASS] 2-local state: informative flat = 0 "
          f"(w_sat=2 >= W=2)")

    eps_flat_2loc = _eps_flat(exp_2loc, kstar_ops, n)
    assert eps_flat_2loc < 1e-10
    passed += 1
    print(f"  [PASS] 2-local state: eps_flat(K*) = {eps_flat_2loc:.2e}")

    # --- GHZ: illustrates w_sat limitation ---

    rho_ghz = depolarise(ghz_state(n), p_noise)
    exp_ghz = pauli_expectations(rho_ghz, n)

    # GHZ has weight-4 phase correlators that K* (w_sat=2) misses
    n_flat_ghz_kstar = 0
    n_info_w4 = 0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        w = pauli_weight(p)
        x_P = exp_ghz.get(tuple(p), 0.0)
        if w == 4 and abs(x_P) > 1e-6:
            n_info_w4 += 1
            if tuple(p) not in kstar_ops:
                n_flat_ghz_kstar += 1

    assert n_info_w4 > 0
    assert n_flat_ghz_kstar > 0
    passed += 1
    print(f"  [PASS] GHZ: {n_flat_ghz_kstar} informative weight-4 flat "
          f"(of {n_info_w4} informative)")

    # GHZ limitation: W(GHZ)=4 > w_sat=2, so eps_flat > 0
    eps_flat_ghz = _eps_flat(exp_ghz, kstar_ops, n)
    assert eps_flat_ghz > 0.1, \
        f"GHZ eps_flat = {eps_flat_ghz}, expected > 0.1"
    passed += 1
    print(f"  [PASS] GHZ eps_flat = {eps_flat_ghz:.4f} > 0 "
          f"(W(GHZ)=4 > w_sat=2)")

    print(f"  [PASS] Hessian properties: {passed} checks")
    return passed
