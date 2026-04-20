"""MLE informative-subspace identifiability (Theorem 1(i))."""

import numpy as np

from common import (
    w_state, product_state_plus, depolarise,
    pauli_tensor, pauli_expectations, all_paulis,
    kstar_operator_indices, random_operator_indices
)
from registry import claims


def _simulate_mle(rho_true, operator_set, n, N_shots=10000, max_iter=500,
                  tol=1e-8):
    """Iterative MLE reconstruction (R*rho*R).  Returns reconstructed rho."""
    d = 2**n
    rng = np.random.default_rng(0)

    measured_ops = []
    measured_freqs = []
    for p_idx in operator_set:
        P = pauli_tensor(p_idx)
        true_exp = np.real(np.trace(P @ rho_true))
        prob_plus = (1 + true_exp) / 2
        n_plus = rng.binomial(N_shots, prob_plus)
        freq = (2 * n_plus / N_shots - 1)
        measured_ops.append(P)
        measured_freqs.append(freq)

    M = len(measured_ops)
    rho = np.eye(d, dtype=complex) / d

    for iteration in range(max_iter):
        R = np.zeros((d, d), dtype=complex)
        for k in range(M):
            P = measured_ops[k]
            f_k = measured_freqs[k]
            pred_k = np.real(np.trace(P @ rho))
            if abs(pred_k) < 1.0:
                R += (f_k - pred_k) * P / M

        eta = 0.5 / (1 + iteration * 0.01)
        update = np.eye(d) + eta * R
        rho_new = update @ rho @ update.conj().T

        eigvals, eigvecs = np.linalg.eigh(rho_new)
        eigvals = np.maximum(eigvals, 0)
        if np.sum(eigvals) > 0:
            eigvals /= np.sum(eigvals)
        else:
            eigvals = np.ones(d) / d
        rho_new = (eigvecs * eigvals) @ eigvecs.conj().T

        diff = np.linalg.norm(rho_new - rho, 'fro')
        rho = rho_new
        if diff < tol:
            break

    return rho


def _fidelity(rho, sigma):
    """F = (tr sqrt(sqrt(rho) sigma sqrt(rho)))^2."""
    eigvals_rho, eigvecs_rho = np.linalg.eigh(rho)
    sqrt_rho = (eigvecs_rho * np.sqrt(np.maximum(eigvals_rho, 0))) \
               @ eigvecs_rho.conj().T
    M = sqrt_rho @ sigma @ sqrt_rho
    eigvals_M = np.linalg.eigvalsh(M)
    return (np.sum(np.sqrt(np.maximum(eigvals_M, 0))))**2


def test_mle_basin():
    """Verify MLE basin properties from Theorem 1(iii)."""
    passed = 0
    n = 4
    p_noise = 0.03
    N_shots = 5000

    # product state: W(rho)=1, perfect structured case
    rho_prod = depolarise(product_state_plus(n), p_noise)
    kstar_ops = kstar_operator_indices(n, K=5)

    rho_kstar = _simulate_mle(rho_prod, kstar_ops, n, N_shots)
    F_kstar_prod = _fidelity(rho_prod, rho_kstar)

    assert F_kstar_prod > 0.90, \
        f"F(K*, product) = {F_kstar_prod:.3f}, expected > 0.90"
    passed += 1
    print(f"  [PASS] Product state F(K*) = {F_kstar_prod:.4f} > 0.90")

    rand_ops = random_operator_indices(n, claims.get("thm:basin", "n4_M"), rng=np.random.default_rng(42))
    rho_rand = _simulate_mle(rho_prod, rand_ops, n, N_shots)
    F_rand_prod = _fidelity(rho_prod, rho_rand)

    assert F_rand_prod > 0.0, \
        f"F(random, product) = {F_rand_prod:.3f}, expected > 0"
    passed += 1
    print(f"  [PASS] Product state F(random) = {F_rand_prod:.4f}")
    print(f"         Delta F = {F_kstar_prod - F_rand_prod:+.4f}")

    # W state
    rho_W = depolarise(w_state(n), p_noise)

    rho_kstar_W = _simulate_mle(rho_W, kstar_ops, n, N_shots)
    F_kstar_W = _fidelity(rho_W, rho_kstar_W)

    rho_rand_W = _simulate_mle(rho_W, rand_ops, n, N_shots)
    F_rand_W = _fidelity(rho_W, rho_rand_W)

    dF_W = F_kstar_W - F_rand_W
    assert dF_W > 0, f"K* should beat random on W state: dF = {dF_W:+.4f}"
    passed += 1
    print(f"  [PASS] W state F(K*) = {F_kstar_W:.4f}, "
          f"F(random) = {F_rand_W:.4f}")
    print(f"         Delta F = {dF_W:+.4f}")

    # --- eps_flat correlation with fidelity advantage ---

    exp_W = pauli_expectations(rho_W, n)

    eps_flat_kstar = 0.0
    eps_flat_rand = 0.0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        x_P = exp_W.get(tuple(p), 0.0)
        if abs(x_P) > 1e-6:
            if tuple(p) not in kstar_ops:
                eps_flat_kstar += x_P**2
            if tuple(p) not in rand_ops:
                eps_flat_rand += x_P**2

    assert eps_flat_kstar < eps_flat_rand, \
        f"eps_flat: K*={eps_flat_kstar:.4f} >= random={eps_flat_rand:.4f}"
    passed += 1
    print(f"  [PASS] eps_flat: K*={eps_flat_kstar:.4f} < "
          f"random={eps_flat_rand:.4f}")

    # K* should win majority across 50 random seeds (statistically powered)
    n_trials = 50
    F_random_all = []
    for seed in range(n_trials):
        rand_seed_ops = random_operator_indices(
            n, claims.get("thm:basin", "n4_M"),
            rng=np.random.default_rng(seed + 100))
        rho_r = _simulate_mle(rho_W, rand_seed_ops, n, N_shots)
        F_r = _fidelity(rho_W, rho_r)
        F_random_all.append(F_r)

    n_wins = sum(1 for F_r in F_random_all if F_kstar_W > F_r)
    F_rand_mean = np.mean(F_random_all)
    F_rand_std = np.std(F_random_all)
    win_rate = n_wins / n_trials

    # Require >= 80% win rate (40/50) for statistical significance
    assert win_rate >= 0.80, \
        f"K* win rate = {win_rate:.0%} ({n_wins}/{n_trials}), expected >= 80%"
    passed += 1
    print(f"  [PASS] K* wins {n_wins}/{n_trials} random seeds "
          f"(win rate = {win_rate:.0%})")
    print(f"         F(random): mean = {F_rand_mean:.4f} "
          f"± {F_rand_std:.4f}")
    print(f"         F(K*) = {F_kstar_W:.4f}, "
          f"mean dF = {F_kstar_W - F_rand_mean:+.4f}")

    print(f"  [PASS] MLE basin: {passed} checks")
    return passed
