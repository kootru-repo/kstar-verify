"""
Shot-Noise-Robust Reconstruction for K* Certification
======================================================
The standard MLE (R*rho*R) collapses under finite shot noise because:
  1. The system is underdetermined (137 < 256 for n=4)
  2. Noisy expectations cause R to have extreme eigenvalues
  3. The iteration amplifies noise instead of converging

This module implements three stabilization mechanisms:
  A. Hedged MLE: regularize toward maximally mixed (prevents rank collapse)
  B. Variance-weighted MLE: weight by inverse shot noise variance
  C. Ratio clipping: b/y ratios clipped to [-5, 5] to prevent divergence
plus damped iteration (alpha * update + (1-alpha) * current state) for
convergence smoothing.

A full Hradil form (reconstruct_robust_mle_full_hradil) is also provided
for analyses where the ratio form's clipping artifacts are problematic
(e.g. GHZ reanalysis in SM Sec. VII D).

Run: PYTHONIOENCODING=utf-8 python robust_mle.py
"""

import numpy as np
from numpy.linalg import norm
from core import (
    compute_kstar, cumulative_count, kstar_measurement_budget,
    select_kstar_paulis, select_random_paulis,
    ghz_state, w_state, random_pure_state, random_mixed_state,
    state_fidelity, project_to_density_matrix,
    depolarizing_channel, readout_error, shot_noise,
)


def reconstruct_robust_mle_full_hradil(expectations, pauli_ops, n_qubits,
                                        n_shots=None, max_iter=1000, tol=1e-7,
                                        hedge=0.05, damping=0.5,
                                        return_iters=False):
    """
    Full Hradil R-operator for the Pauli binomial likelihood.

    R = sum_i w_i * [(1-b_i*y_i)/(1-y_i^2) * I + (b_i-y_i)/(1-y_i^2) * P_i]

    Both coefficients are regular at y=0 and R=I at the fixed point b=y.
    This form avoids the clipping artifacts of the ratio form, but its
    identity contribution slows convergence in underdetermined regimes.
    Used for the GHZ reanalysis (SM Sec. VII D) where it recovers F=0.888
    vs the ratio form's F=0.496.

    Args:
        expectations: measured Pauli expectation values
        pauli_ops: list of Pauli operator matrices
        n_qubits: number of qubits
        n_shots: shot count (None = infinite, no variance weighting)
        max_iter: maximum iterations
        tol: convergence tolerance (Frobenius norm change)
        hedge: regularization strength toward maximally mixed (0.01-0.1)
        damping: initial damping factor (0 = no update, 1 = full R*rho*R)
    """
    dim = 2 ** n_qubits
    m = len(expectations)
    b = np.array(expectations, dtype=float)

    # Compute per-measurement weights
    if n_shots is not None and n_shots > 0:
        # Variance of Pauli expectation from binomial: var = (1 - <P>^2) / N
        # Clip expectations to avoid division by zero
        b_clipped = np.clip(b, -0.999, 0.999)
        variances = (1.0 - b_clipped ** 2) / n_shots
        # Weight = 1/variance, normalized
        weights = 1.0 / np.maximum(variances, 1e-10)
        weights /= weights.sum()
        weights *= m  # Scale so sum = m (same as unweighted)
    else:
        weights = np.ones(m)

    # Initialize with maximally mixed state
    rho = np.eye(dim, dtype=complex) / dim
    identity = np.eye(dim, dtype=complex) / dim

    # Precompute
    P_ops = [P.copy() for P in pauli_ops]

    alpha = damping  # Start with cautious updates

    for iteration in range(max_iter):
        rho_old = rho.copy()

        # Predicted expectations
        predicted = np.array([np.trace(P @ rho).real for P in P_ops])

        # Build weighted R operator with hedge (full Hradil form)
        R = hedge * np.eye(dim, dtype=complex)
        for i in range(m):
            y_i = predicted[i]
            denom = max(1.0 - y_i**2, 1e-12)
            alpha_i = (1.0 - b[i] * y_i) / denom
            beta_i = (b[i] - y_i) / denom
            R += weights[i] * (alpha_i * np.eye(dim, dtype=complex)
                               + beta_i * P_ops[i])
        R /= (m + hedge * dim)

        # Damped R*rho*R iteration
        rho_update = R @ rho @ R.conj().T
        # Ensure Hermitian and normalize
        rho_update = (rho_update + rho_update.conj().T) / 2
        tr = np.trace(rho_update).real
        if tr > 1e-15:
            rho_update /= tr

        # Damped update: mix new and old
        rho = alpha * rho_update + (1 - alpha) * rho

        # Re-normalize
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > 1e-15:
            rho /= tr

        # Convergence check
        diff = norm(rho - rho_old, 'fro')
        if diff < tol:
            break

        # Adaptive damping: increase alpha as we converge
        if iteration > 20 and diff < 0.01:
            alpha = min(alpha * 1.05, 0.95)

    # Final projection to ensure valid density matrix
    rho = project_to_density_matrix(rho)
    if return_iters:
        return rho, iteration + 1
    return rho


def reconstruct_robust_mle(expectations, pauli_ops, n_qubits,
                           n_shots=None, max_iter=1000, tol=1e-7,
                           hedge=0.05, damping=0.5, return_iters=False):
    """
    Hedged, variance-weighted, damped MLE reconstruction (ratio form).

    R = sum_i w_i * clip(b_i/y_i, -5, 5) * P_i

    This is the primary reconstruction algorithm used for all hardware
    results reported in the paper (Table I, Table III, etc.). The ratio
    form with clipping is a standard iterative MLE approach (Hradil 1997,
    Paris 2004).
    """
    dim = 2 ** n_qubits
    m = len(expectations)
    b = np.array(expectations, dtype=float)

    if n_shots is not None and n_shots > 0:
        b_clipped = np.clip(b, -0.999, 0.999)
        variances = (1.0 - b_clipped ** 2) / n_shots
        weights = 1.0 / np.maximum(variances, 1e-10)
        weights /= weights.sum()
        weights *= m
    else:
        weights = np.ones(m)

    rho = np.eye(dim, dtype=complex) / dim
    P_ops = [P.copy() for P in pauli_ops]
    alpha = damping

    for iteration in range(max_iter):
        rho_old = rho.copy()
        predicted = np.array([np.trace(P @ rho).real for P in P_ops])

        R = hedge * np.eye(dim, dtype=complex)
        for i in range(m):
            pred_i = predicted[i]
            if abs(pred_i) < 1e-12:
                pred_i = np.sign(pred_i + 1e-20) * 1e-12
            ratio = b[i] / pred_i
            ratio = np.clip(ratio, -5.0, 5.0)
            R += weights[i] * ratio * P_ops[i]
        R /= (m + hedge * dim)

        rho_update = R @ rho @ R.conj().T
        rho_update = (rho_update + rho_update.conj().T) / 2
        tr = np.trace(rho_update).real
        if tr > 1e-15:
            rho_update /= tr

        rho = alpha * rho_update + (1 - alpha) * rho
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > 1e-15:
            rho /= tr

        diff = norm(rho - rho_old, 'fro')
        if diff < tol:
            break
        if iteration > 20 and diff < 0.01:
            alpha = min(alpha * 1.05, 0.95)

    rho = project_to_density_matrix(rho)
    if return_iters:
        return rho, iteration + 1
    return rho


def reconstruct_hedged_lstsq(expectations, pauli_ops, n_qubits,
                              n_shots=None, hedge=0.1):
    """
    Tikhonov-regularized least squares with variance weighting.

    Solves: min ||W(A*rho_vec - b)||^2 + lambda*||rho_vec - rho_mm||^2
    where W = diag(1/sigma_i), rho_mm = I/d (maximally mixed).

    This is simpler and more stable than iterative MLE for noisy data.
    """
    dim = 2 ** n_qubits
    d2 = dim * dim
    m = len(expectations)
    b = np.array(expectations, dtype=float)

    # Build measurement matrix
    A = np.zeros((m, d2), dtype=complex)
    for i, P in enumerate(pauli_ops):
        A[i, :] = P.flatten()

    # Variance weighting
    if n_shots is not None and n_shots > 0:
        b_clipped = np.clip(b, -0.999, 0.999)
        sigma = np.sqrt((1.0 - b_clipped ** 2) / n_shots)
        sigma = np.maximum(sigma, 1e-6)
        W = np.diag(1.0 / sigma)
        A_w = W @ A
        b_w = W @ b
    else:
        A_w = A
        b_w = b

    # Tikhonov regularization toward maximally mixed
    rho_mm = np.eye(dim, dtype=complex).flatten() / dim
    lam = hedge * m  # Regularization strength scales with data size

    # Augmented system: [A_w; sqrt(lam)*I] @ x = [b_w; sqrt(lam)*rho_mm]
    A_aug = np.vstack([A_w, np.sqrt(lam) * np.eye(d2, dtype=complex)])
    b_aug = np.concatenate([b_w, np.sqrt(lam) * rho_mm])

    rho_vec, _, _, _ = np.linalg.lstsq(A_aug, b_aug, rcond=None)
    return project_to_density_matrix(rho_vec.reshape(dim, dim))


# ====================================================================
# Test harness
# ====================================================================

def test_reconstruction(n_qubits, rho_true, state_name, pauli_ops,
                        selection_name, noise_model=None, n_shots=None):
    """Test all reconstruction methods on one configuration."""
    dim = 2 ** n_qubits

    # Apply state-level noise
    rho_noisy = rho_true.copy()
    if noise_model and 'depolarizing' in noise_model:
        p = noise_model['depolarizing']
        rho_noisy = (1 - p) * rho_noisy + p * np.eye(dim) / dim

    # Measure
    expectations = []
    for P in pauli_ops:
        exp = np.trace(rho_noisy @ P).real
        if noise_model and 'readout_error' in noise_model:
            exp = readout_error(exp, noise_model['readout_error'])
        if n_shots is not None:
            exp = shot_noise(exp, n_shots)
        expectations.append(exp)

    results = {}

    # Standard MLE
    from prove_selection_advantage import reconstruct_mle
    rho_std = reconstruct_mle(expectations, pauli_ops, n_qubits)
    results['standard_mle'] = state_fidelity(rho_true, rho_std)

    # Robust MLE
    rho_rob = reconstruct_robust_mle(expectations, pauli_ops, n_qubits,
                                      n_shots=n_shots)
    results['robust_mle'] = state_fidelity(rho_true, rho_rob)

    # Hedged lstsq
    rho_hlsq = reconstruct_hedged_lstsq(expectations, pauli_ops, n_qubits,
                                         n_shots=n_shots)
    results['hedged_lstsq'] = state_fidelity(rho_true, rho_hlsq)

    return results


def main():
    np.random.seed(42)

    print("=" * 85)
    print("  SHOT-NOISE-ROBUST RECONSTRUCTION — ALGORITHM COMPARISON")
    print("  Testing: standard MLE vs robust MLE vs hedged lstsq")
    print("=" * 85)

    # ================================================================
    # Test 1: GHZ(n=3) with K* selection under various shot counts
    # ================================================================
    print(f"\n{'='*85}")
    print(f"  TEST 1: GHZ(n=3), K*-selected Paulis, varying shot count")
    print(f"{'='*85}")

    n = 3
    psi = ghz_state(n)
    rho = np.outer(psi, psi.conj())
    kops, _, _ = select_kstar_paulis(n)

    print(f"\n  {'Shots':>10s}  {'Std MLE':>10s}  {'Robust MLE':>10s}  {'Hedged LS':>10s}  {'Best':>12s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*12}")

    for shots in [None, 100000, 10000, 4000, 1000, 500, 100]:
        # Average over multiple noise draws
        fids = {'standard_mle': [], 'robust_mle': [], 'hedged_lstsq': []}
        n_draws = 1 if shots is None else 10

        for _ in range(n_draws):
            r = test_reconstruction(n, rho, 'GHZ', kops, 'K*', n_shots=shots)
            for k, v in r.items():
                fids[k].append(v)

        means = {k: np.mean(v) for k, v in fids.items()}
        best = max(means, key=means.get)
        shots_str = 'inf' if shots is None else str(shots)
        print(f"  {shots_str:>10s}  {means['standard_mle']:>10.4f}  "
              f"{means['robust_mle']:>10.4f}  {means['hedged_lstsq']:>10.4f}  "
              f"{best:>12s}")

    # ================================================================
    # Test 2: GHZ(n=4) with K* selection under various shot counts
    # ================================================================
    print(f"\n{'='*85}")
    print(f"  TEST 2: GHZ(n=4), K*-selected Paulis, varying shot count")
    print(f"{'='*85}")

    n = 4
    psi = ghz_state(n)
    rho = np.outer(psi, psi.conj())
    kops, _, _ = select_kstar_paulis(n)

    print(f"\n  {'Shots':>10s}  {'Std MLE':>10s}  {'Robust MLE':>10s}  {'Hedged LS':>10s}  {'Best':>12s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*12}")

    for shots in [None, 100000, 10000, 4000, 1000, 500, 100]:
        fids = {'standard_mle': [], 'robust_mle': [], 'hedged_lstsq': []}
        n_draws = 1 if shots is None else 8

        for _ in range(n_draws):
            r = test_reconstruction(n, rho, 'GHZ', kops, 'K*', n_shots=shots)
            for k, v in r.items():
                fids[k].append(v)

        means = {k: np.mean(v) for k, v in fids.items()}
        best = max(means, key=means.get)
        shots_str = 'inf' if shots is None else str(shots)
        print(f"  {shots_str:>10s}  {means['standard_mle']:>10.4f}  "
              f"{means['robust_mle']:>10.4f}  {means['hedged_lstsq']:>10.4f}  "
              f"{best:>12s}")

    # ================================================================
    # Test 3: GHZ(n=4) + combined noise (the failure case)
    # ================================================================
    print(f"\n{'='*85}")
    print(f"  TEST 3: GHZ(n=4), K*-selected, combined noise scenarios")
    print(f"{'='*85}")

    n = 4
    psi = ghz_state(n)
    rho = np.outer(psi, psi.conj())
    kops, _, _ = select_kstar_paulis(n)

    noise_configs = [
        ('noiseless',         None,                                        None),
        ('depol_5pct',        {'depolarizing': 0.05},                      None),
        ('depol_5pct+4k',     {'depolarizing': 0.05},                      4000),
        ('depol_5pct+1k',     {'depolarizing': 0.05},                      1000),
        ('rdout_2pct+4k',     {'readout_error': 0.02},                     4000),
        ('rdout_2pct+1k',     {'readout_error': 0.02},                     1000),
        ('mild_comb+4k',      {'depolarizing': 0.01, 'readout_error': 0.02}, 4000),
        ('mild_comb+1k',      {'depolarizing': 0.01, 'readout_error': 0.02}, 1000),
        ('mod_comb+1k',       {'depolarizing': 0.05, 'readout_error': 0.05}, 1000),
        ('harsh_comb+100',    {'depolarizing': 0.10, 'readout_error': 0.05}, 100),
    ]

    print(f"\n  {'Noise':20s}  {'Std MLE':>10s}  {'Robust MLE':>10s}  {'Hedged LS':>10s}  {'Best':>12s}")
    print(f"  {'-'*20}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*12}")

    for label, noise, shots in noise_configs:
        fids = {'standard_mle': [], 'robust_mle': [], 'hedged_lstsq': []}
        n_draws = 1 if shots is None else 8

        for _ in range(n_draws):
            r = test_reconstruction(n, rho, 'GHZ', kops, 'K*',
                                    noise_model=noise, n_shots=shots)
            for k, v in r.items():
                fids[k].append(v)

        means = {k: np.mean(v) for k, v in fids.items()}
        best = max(means, key=means.get)
        print(f"  {label:20s}  {means['standard_mle']:>10.4f}  "
              f"{means['robust_mle']:>10.4f}  {means['hedged_lstsq']:>10.4f}  "
              f"{best:>12s}")

    # ================================================================
    # Test 4: K* vs Random with ROBUST reconstruction
    # The key question: does K* advantage survive noise with better recon?
    # ================================================================
    print(f"\n{'='*85}")
    print(f"  TEST 4: K* vs RANDOM with robust MLE — does advantage survive?")
    print(f"{'='*85}")

    for n in [3, 4]:
        psi = ghz_state(n)
        rho = np.outer(psi, psi.conj())
        budget = kstar_measurement_budget(n)
        n_meas = budget['n_modes']

        kops, _, _ = select_kstar_paulis(n)

        print(f"\n  --- GHZ(n={n}), budget={n_meas} ---")
        print(f"  {'Noise':20s}  {'K* Rob':>8s}  {'Rand Rob':>10s}  {'Adv':>8s}  {'Win?':>5s}")
        print(f"  {'-'*20}  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*5}")

        test_configs = [
            ('noiseless',       None,                                        None),
            ('shots_4000',      None,                                        4000),
            ('shots_1000',      None,                                        1000),
            ('depol_5pct+4k',   {'depolarizing': 0.05},                      4000),
            ('depol_5pct+1k',   {'depolarizing': 0.05},                      1000),
            ('mild_comb+4k',    {'depolarizing': 0.01, 'readout_error': 0.02}, 4000),
            ('mild_comb+1k',    {'depolarizing': 0.01, 'readout_error': 0.02}, 1000),
            ('mod_comb+1k',     {'depolarizing': 0.05, 'readout_error': 0.05}, 1000),
        ]

        for label, noise, shots in test_configs:
            dim = 2 ** n

            # K* with robust MLE (multiple draws)
            n_draws = 1 if shots is None else 8
            k_fids = []
            for _ in range(n_draws):
                rho_noisy = rho.copy()
                if noise and 'depolarizing' in noise:
                    p = noise['depolarizing']
                    rho_noisy = (1 - p) * rho_noisy + p * np.eye(dim) / dim
                kexp = []
                for P in kops:
                    exp = np.trace(rho_noisy @ P).real
                    if noise and 'readout_error' in noise:
                        exp = readout_error(exp, noise['readout_error'])
                    if shots is not None:
                        exp = shot_noise(exp, shots)
                    kexp.append(exp)
                rho_k = reconstruct_robust_mle(kexp, kops, n, n_shots=shots)
                k_fids.append(state_fidelity(rho, rho_k))

            # Random with robust MLE (multiple seeds x draws)
            n_seeds = 10 if n <= 3 else 6
            r_fids = []
            for seed in range(n_seeds):
                rops, _, _ = select_random_paulis(n, n_meas, seed=seed)
                for _ in range(max(1, n_draws // 2)):
                    rho_noisy = rho.copy()
                    if noise and 'depolarizing' in noise:
                        p = noise['depolarizing']
                        rho_noisy = (1 - p) * rho_noisy + p * np.eye(dim) / dim
                    rexp = []
                    for P in rops:
                        exp = np.trace(rho_noisy @ P).real
                        if noise and 'readout_error' in noise:
                            exp = readout_error(exp, noise['readout_error'])
                        if shots is not None:
                            exp = shot_noise(exp, shots)
                        rexp.append(exp)
                    rho_r = reconstruct_robust_mle(rexp, rops, n, n_shots=shots)
                    r_fids.append(state_fidelity(rho, rho_r))

            kf = np.mean(k_fids)
            rf_mean = np.mean(r_fids)
            rf_std = np.std(r_fids)
            adv = kf - rf_mean
            wins = kf > rf_mean + rf_std

            win_str = "YES" if wins else "no"
            print(f"  {label:20s}  {kf:>8.4f}  "
                  f"{rf_mean:>7.4f}+/-{rf_std:.3f}  "
                  f"{adv:>+8.4f}  {win_str:>5s}")

    # ================================================================
    # Test 5: Random pure states — control (should NOT show K* advantage)
    # ================================================================
    print(f"\n{'='*85}")
    print(f"  TEST 5: Random pure states (n=3) — control group")
    print(f"{'='*85}")

    n = 3
    budget = kstar_measurement_budget(n)
    n_meas = budget['n_modes']
    kops, _, _ = select_kstar_paulis(n)

    print(f"\n  {'State':15s}  {'Noise':15s}  {'K* Rob':>8s}  {'Rand Rob':>8s}  {'Adv':>8s}")
    print(f"  {'-'*15}  {'-'*15}  {'-'*8}  {'-'*8}  {'-'*8}")

    for i in range(5):
        psi = random_pure_state(n)
        rho_test = np.outer(psi, psi.conj())

        for label, noise, shots in [('noiseless', None, None),
                                     ('shots_1000', None, 1000)]:
            dim = 2 ** n
            n_draws = 1 if shots is None else 5

            k_fids = []
            for _ in range(n_draws):
                rho_noisy = rho_test.copy()
                kexp = []
                for P in kops:
                    exp = np.trace(rho_noisy @ P).real
                    if shots: exp = shot_noise(exp, shots)
                    kexp.append(exp)
                rho_k = reconstruct_robust_mle(kexp, kops, n, n_shots=shots)
                k_fids.append(state_fidelity(rho_test, rho_k))

            r_fids = []
            for seed in range(8):
                rops, _, _ = select_random_paulis(n, n_meas, seed=seed)
                for _ in range(max(1, n_draws // 2)):
                    rexp = []
                    for P in rops:
                        exp = np.trace(rho_noisy @ P).real
                        if shots: exp = shot_noise(exp, shots)
                        rexp.append(exp)
                    rho_r = reconstruct_robust_mle(rexp, rops, n, n_shots=shots)
                    r_fids.append(state_fidelity(rho_test, rho_r))

            kf = np.mean(k_fids)
            rf = np.mean(r_fids)
            print(f"  rand_pure_{i:3d}  {label:15s}  {kf:>8.4f}  {rf:>8.4f}  {kf-rf:>+8.4f}")

    print(f"\n{'='*85}")
    print(f"  DONE")
    print(f"{'='*85}")


if __name__ == '__main__':
    main()
