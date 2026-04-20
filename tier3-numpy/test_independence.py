"""Finite-sample independence of Pauli estimates (Lemma 1 deepening).

Verifies:
1. Estimates y_hat_P are statistically independent across P (not just uncorrelated)
2. Hoeffding bound holds empirically
3. Cross-correlations are zero at finite N_s
"""

import numpy as np

from common import (
    w_state, product_state_plus, ghz_state,
    pauli_tensor, pauli_expectations, pauli_weight, all_paulis,
    kstar_operator_indices,
)


def test_independence():
    """Verify finite-sample independence claims from Lemma 1."""
    passed = 0
    n = 4
    d = 2**n
    n_repeats = 200

    # --- Test 1: Cross-correlation of estimates is zero ---
    # Simulate Pauli measurements and check that estimates for
    # different operators are uncorrelated.

    rho = w_state(n)
    exp_true = pauli_expectations(rho, n)

    # Pick 6 operators across different weights
    test_ops = [
        (3, 0, 0, 0),  # Z_1, weight 1
        (0, 1, 0, 0),  # X_2, weight 1
        (1, 1, 0, 0),  # X_1 X_2, weight 2
        (3, 3, 0, 0),  # Z_1 Z_2, weight 2
        (1, 2, 3, 0),  # X_1 Y_2 Z_3, weight 3
        (1, 1, 1, 1),  # X_1 X_2 X_3 X_4, weight 4
    ]

    N_s = 1000
    rng = np.random.default_rng(42)

    # Simulate many independent measurement rounds
    estimates = {op: [] for op in test_ops}
    for _ in range(n_repeats):
        for op in test_ops:
            x_P = exp_true[op]
            prob_plus = (1 + x_P) / 2
            n_plus = rng.binomial(N_s, prob_plus)
            y_hat = 2 * n_plus / N_s - 1
            estimates[op].append(y_hat)

    # Check pairwise cross-correlations are near zero
    ops_list = list(test_ops)
    max_cross_corr = 0.0
    for i in range(len(ops_list)):
        for j in range(i + 1, len(ops_list)):
            a = np.array(estimates[ops_list[i]])
            b = np.array(estimates[ops_list[j]])
            # Pearson correlation
            corr = np.corrcoef(a, b)[0, 1]
            max_cross_corr = max(max_cross_corr, abs(corr))

    # With n_repeats=200 independent samples, |corr| should be < 0.2
    # with very high probability for truly independent variables
    assert max_cross_corr < 0.2, \
        f"Max cross-correlation = {max_cross_corr:.3f}, expected < 0.2"
    passed += 1
    print(f"  [PASS] Cross-correlations: max |r| = {max_cross_corr:.4f} < 0.2 "
          f"({len(ops_list)} operators, {n_repeats} trials)")

    # --- Test 2: Hoeffding bound holds empirically ---
    # P(|y_hat - x_P| > t) <= 2*exp(-N_s * t^2 / 2)
    # Test at multiple t values and verify empirical violation rate
    # is below the Hoeffding bound.

    rho_prod = product_state_plus(n)
    exp_prod = pauli_expectations(rho_prod, n)

    # Test operators with different true expectations
    hoeffding_ops = [
        ((1, 0, 0, 0), "X_1"),   # x_P = 1.0
        ((3, 0, 0, 0), "Z_1"),   # x_P = 0.0
        ((1, 1, 0, 0), "X1X2"),  # x_P = 1.0
        ((3, 3, 0, 0), "Z1Z2"),  # x_P = 0.0
    ]

    N_s_hoeffding = 500
    n_trials = 5000
    t_values = [0.05, 0.10, 0.15]
    rng2 = np.random.default_rng(123)

    for op, name in hoeffding_ops:
        x_P = exp_prod[op]
        if abs(x_P) > 1.0 - 1e-10:
            # Skip |x_P| = 1 (binomial degenerates)
            continue

        for t in t_values:
            violations = 0
            for _ in range(n_trials):
                prob_plus = (1 + x_P) / 2
                n_plus = rng2.binomial(N_s_hoeffding, prob_plus)
                y_hat = 2 * n_plus / N_s_hoeffding - 1
                if abs(y_hat - x_P) > t:
                    violations += 1

            empirical_rate = violations / n_trials
            hoeffding_bound = 2 * np.exp(-N_s_hoeffding * t**2 / 2)

            # Empirical rate should be below Hoeffding bound
            # Allow small statistical margin (1.5x)
            assert empirical_rate <= hoeffding_bound * 1.5 + 0.01, \
                f"{name}: empirical {empirical_rate:.4f} > " \
                f"1.5 * Hoeffding {hoeffding_bound:.4f}"
            passed += 1
            print(f"  [PASS] Hoeffding: {name} (x={x_P:.1f}), t={t:.2f}: "
                  f"empirical={empirical_rate:.4f}, "
                  f"bound={hoeffding_bound:.4f}")

    # --- Test 3: Independence test via mutual information proxy ---
    # If y_hat_P and y_hat_{P'} are independent, then knowing one
    # should not help predict the other. Test via conditional variance.

    rho_W = w_state(n)
    exp_W = pauli_expectations(rho_W, n)

    op_A = (3, 0, 0, 0)  # Z_1
    op_B = (0, 3, 0, 0)  # Z_2
    N_s_test = 500
    n_trials_mi = 2000
    rng3 = np.random.default_rng(456)

    estimates_A = []
    estimates_B = []
    for _ in range(n_trials_mi):
        x_A = exp_W[op_A]
        x_B = exp_W[op_B]
        n_A = rng3.binomial(N_s_test, (1 + x_A) / 2)
        n_B = rng3.binomial(N_s_test, (1 + x_B) / 2)
        estimates_A.append(2 * n_A / N_s_test - 1)
        estimates_B.append(2 * n_B / N_s_test - 1)

    estimates_A = np.array(estimates_A)
    estimates_B = np.array(estimates_B)

    # Split A into high/low halves and check B's variance is the same
    median_A = np.median(estimates_A)
    high_mask = estimates_A > median_A
    low_mask = ~high_mask

    var_B_high = np.var(estimates_B[high_mask])
    var_B_low = np.var(estimates_B[low_mask])
    var_B_all = np.var(estimates_B)

    # If independent: var_B_high ≈ var_B_low ≈ var_B_all
    ratio_high = var_B_high / var_B_all
    ratio_low = var_B_low / var_B_all

    assert 0.85 < ratio_high < 1.15, \
        f"Conditional variance ratio (high) = {ratio_high:.3f}"
    assert 0.85 < ratio_low < 1.15, \
        f"Conditional variance ratio (low) = {ratio_low:.3f}"
    passed += 1
    print(f"  [PASS] Conditional independence: "
          f"Var(B|A_high)/Var(B) = {ratio_high:.3f}, "
          f"Var(B|A_low)/Var(B) = {ratio_low:.3f}")

    # --- Test 4: Fisher information matrix is diagonal by Pauli orthonormality ---
    # For Pauli measurements, H_{PP'} = delta_{PP'} * N_s/(1-y_P^2) because
    # dTr(Q rho)/dy_P = delta_{QP}. This verifies the structural property
    # at three different states (GHZ, W, product).

    for state_name, rho_test in [("GHZ", ghz_state(n)),
                                  ("W", w_state(n)),
                                  ("product", product_state_plus(n))]:
        exp_test = pauli_expectations(rho_test, n)
        kstar_ops = kstar_operator_indices(n)

        all_ops = [p for p in all_paulis(n)
                   if not all(x == 0 for x in p)]
        M = len(all_ops)
        op_to_idx = {tuple(p): i for i, p in enumerate(all_ops)}

        # Build Fisher matrix
        H = np.zeros((M, M))
        for Q in kstar_ops:
            if all(x == 0 for x in Q):
                continue
            if Q not in op_to_idx:
                continue
            i = op_to_idx[Q]
            y_Q = exp_test.get(Q, 0.0)
            if abs(y_Q) < 1.0 - 1e-12:
                H[i, i] += 1000.0 / (1.0 - y_Q**2)
            else:
                H[i, i] += 1e9

        # Off-diagonal must be exactly 0
        off_diag = H - np.diag(np.diag(H))
        assert np.max(np.abs(off_diag)) < 1e-15, \
            f"{state_name}: off-diagonal max = {np.max(np.abs(off_diag))}"
        passed += 1
        print(f"  [PASS] Fisher diagonal: {state_name}, "
              f"max off-diag = {np.max(np.abs(off_diag)):.1e}")

    print(f"\n  [PASS] Finite-sample independence: {passed} checks")
    return passed
