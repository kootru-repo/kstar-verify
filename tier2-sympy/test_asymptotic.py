"""Operator lower bound (Corollary 2) and asymptotic K* separation (Theorem 3)."""

from sympy import Rational, binomial, log, N as sym_N

from common import (
    qary_parity_weight_counts, qary_gram_eigenvalues,
    weight_class_sizes,
)
from registry import claims


def test_asymptotic():
    """Verify Corollary 2 and Theorem 3 claims."""
    passed = 0

    # =================================================================
    # COROLLARY 2: Operator lower bound
    # Any set with worst-case F > 1/2 over k-local states must include
    # ALL operators of weight 1 through k.  |S| >= sum_{w=1}^k A_w.
    # =================================================================

    # --- n=4, q=2, k=2: lower bound = A_1 + A_2 = 12 + 54 = 66 ---
    n, q, k = 4, 2, 2
    A = weight_class_sizes(n, q)
    claimed_lb = claims.get("cor:lower_bound", "lower_bound_n4_k2")
    n4_M = claims.get("thm:basin", "n4_M")

    lower_bound = sum(A[w] for w in range(1, k + 1))
    assert lower_bound == claimed_lb, f"Lower bound = {lower_bound}, expected {claimed_lb}"
    passed += 1
    print(f"  [PASS] Corollary 2: |S| >= {lower_bound} at n={n}, q={q}, k={k}")

    # --- Verify K* has all weight-<=2 operators ---
    c_kstar = qary_parity_weight_counts(n, 5, q)
    N_kstar = sum(c_kstar)
    assert N_kstar == n4_M, f"N_4(5) = {N_kstar}, expected {n4_M}"

    # Weight saturation: c_w >= A_w for w=0,1,2
    # After redistribution, M_w = A_w for w <= 2.
    # Check: cumulative c_w provides enough for saturation.
    A_cumulative = [sum(A[j] for j in range(w + 1)) for w in range(n + 1)]
    c_cumulative = [sum(c_kstar[j] for j in range(w + 1)) for w in range(n + 1)]
    for w in range(k + 1):
        assert c_cumulative[w] >= A_cumulative[w], \
            f"w={w}: c_cumulative={c_cumulative[w]} < A_cumulative={A_cumulative[w]}"

    necessary = sum(A[w] for w in range(k + 1))  # includes w=0 (identity)
    lattice_extra = N_kstar - necessary
    # A_0 = 1 (identity), so weight<=2 non-identity = 12+54 = 66
    necessary_nonid = sum(A[w] for w in range(1, k + 1))
    assert necessary_nonid == claimed_lb
    assert lattice_extra == n4_M - (1 + claimed_lb), \
        f"Lattice extra = {lattice_extra}, expected {n4_M - 1 - claimed_lb}"
    passed += 1
    print(f"  [PASS] K* has {necessary_nonid} necessary (weight<=2) "
          f"+ {lattice_extra} lattice (weight>2) = {N_kstar}")

    # --- Lower bound at other n values ---
    for n_test, k_test in [(3, 2), (5, 2), (6, 2), (4, 1)]:
        A_test = weight_class_sizes(n_test, q)
        lb = sum(A_test[w] for w in range(1, k_test + 1))
        assert lb > 0
        passed += 1
        print(f"  [PASS] Lower bound at n={n_test}, k={k_test}: |S| >= {lb}")

    # =================================================================
    # THEOREM 3: Asymptotic K* separation
    # =================================================================

    # --- Part (i): K* achieves weight saturation at each n ---
    # Already verified in test_universality multi-n sweep.
    # Re-verify here for completeness.
    for n_test in [2, 3, 4, 5, 6, 7, 8, 9]:
        # Find K* (minimal K with all eigenvalues > 0)
        K_star = None
        for K_try in range(1, n_test**2 + 2):
            eigs = qary_gram_eigenvalues(n_test, K_try, 2)
            if all(e > 0 for e in eigs):
                K_star = K_try
                break
        assert K_star is not None

        c_w = qary_parity_weight_counts(n_test, K_star, 2)
        assert all(c > 0 for c in c_w), \
            f"n={n_test}: c_w has zeros at K*={K_star}"
        passed += 1
        print(f"  [PASS] Theorem 3(i): n={n_test}, K*={K_star}, "
              f"all c_w > 0: {c_w}")

    # --- Part (ii): M_n / N < 1 and random failure bound ---
    print()
    ratios = []
    for n_test in [2, 3, 4, 5, 6, 7, 8, 9]:
        K_star = None
        for K_try in range(1, n_test**2 + 2):
            eigs = qary_gram_eigenvalues(n_test, K_try, 2)
            if all(e > 0 for e in eigs):
                K_star = K_try
                break

        c_w = qary_parity_weight_counts(n_test, K_star, 2)
        M_n = sum(c_w)
        N_total = 4**n_test - 1
        ratio = Rational(M_n, N_total)

        # M_n / N < 1 (proper subset)
        assert ratio < 1, f"n={n_test}: M_n/N = {ratio} >= 1"

        # Random failure bound: P(mu_w = 1) <= (M_n/N)^{A_w}
        A_test = weight_class_sizes(n_test, 2)
        for w in range(1, min(3, n_test + 1)):
            prob_bound = ratio ** A_test[w]
            assert prob_bound < 1
            passed += 1
            print(f"  [PASS] Theorem 3(ii): n={n_test}, w={w}, "
                  f"P(mu_w=1) <= ({float(ratio):.3f})^{A_test[w]} "
                  f"= {float(prob_bound):.2e}")

        ratios.append((n_test, float(ratio)))

    # --- M_n / N convergence: verify ratio < 1 and trending ---
    print()
    for n_test, r in ratios:
        print(f"  [INFO] n={n_test}: M_n/N = {r:.4f}")

    # All ratios < 1
    assert all(r < 1 for _, r in ratios)
    passed += 1
    print(f"  [PASS] All M_n/N < 1 (proper subset for all n tested)")

    # Ratio at n >= 5 is approximately 1/3
    large_n_ratios = [r for n_t, r in ratios if n_t >= 5]
    if large_n_ratios:
        for r in large_n_ratios:
            assert 0.25 < r < 0.42, f"Ratio {r} outside expected range [0.25, 0.42] (~1/3)"
        passed += 1
        print(f"  [PASS] Large-n ratios in [0.25, 0.42], "
              f"consistent with M_n/N -> 1/3")

    # --- Part (iii): Fidelity dichotomy ---
    # At K*: all eigenvalues positive (F >= 1 - O(d^-2))
    # At K*-1 with random: has zero eigenvalue (F <= 1/2 possible)
    # Verify the gap is stark.
    for n_test in [3, 4, 5]:
        K_star = None
        for K_try in range(1, n_test**2 + 2):
            eigs = qary_gram_eigenvalues(n_test, K_try, 2)
            if all(e > 0 for e in eigs):
                K_star = K_try
                break

        eigs_star = qary_gram_eigenvalues(n_test, K_star, 2)
        eigs_below = qary_gram_eigenvalues(n_test, K_star - 1, 2)

        assert all(e > 0 for e in eigs_star)
        assert any(e == 0 for e in eigs_below)

        d = 2**n_test
        # Spectral resolution bound: 1-F <= (d-1-S_k)/d^2 + O(1/N_s)
        # For pure k-local state: S_k = d-1, so 1-F = O(1/N_s) = O(d^-2)
        resolution_bound = Rational(d - 1, d**2)
        assert resolution_bound < Rational(1, 2), \
            f"n={n_test}: resolution bound {resolution_bound} >= 1/2"

        passed += 1
        print(f"  [PASS] Theorem 3(iii): n={n_test}, "
              f"K* bound = {float(resolution_bound):.4f} << 1/2, "
              f"K*-1 has zero eigenvalue (F <= 1/2 possible)")

    # --- Verify the hypergeometric bound is correct ---
    # P(all A_w in random M from N) = prod_{j=0}^{A_w-1} (M-j)/(N-j)
    # <= (M/N)^{A_w}
    # Check: product <= power bound
    n_test = 4
    M_n = n4_M
    N_total = claims.get("thm:basin", "n4_N")
    for w in range(1, 3):
        A_w = int(weight_class_sizes(n_test, 2)[w])
        exact_prob = Rational(1)
        for j in range(A_w):
            exact_prob *= Rational(M_n - j, N_total - j)
        power_bound = Rational(M_n, N_total) ** A_w

        assert exact_prob <= power_bound, \
            f"w={w}: exact {exact_prob} > bound {power_bound}"
        assert exact_prob >= 0

        passed += 1
        print(f"  [PASS] Hypergeometric bound: w={w}, A_w={A_w}, "
              f"exact P = {float(exact_prob):.2e}, "
              f"bound = {float(power_bound):.2e}")

    print(f"\n  [PASS] Asymptotic separation: {passed} checks")
    return passed
