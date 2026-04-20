"""
Domain: Number Theory
=====================
Divisor sums, lattice point counts, Jacobi theta functions.
Portable utilities.
All computations exact over ZZ/QQ.
"""

# Framework (register_test, certify, QQ, ZZ, etc.) loaded by run_all.sage

from itertools import product as cartesian

# ---------------------------------------------------------------------------
# Divisor sums
# ---------------------------------------------------------------------------

def sigma_tilde(k):
    """Sum of divisors of k not divisible by 4. Exact integer."""
    if k == 0:
        return ZZ(0)
    return sum(d for d in divisors(k) if d % 4 != 0)


def r4_jacobi(k):
    """r_4(k) = 8 * sigma_tilde(k) for k >= 1; r_4(0) = 1.
    Jacobi's four-square theorem."""
    if k == 0:
        return ZZ(1)
    return 8 * sigma_tilde(k)


def r4_direct(k):
    """Direct brute-force count of vectors n in Z^4 with |n|^2 = k."""
    bound = isqrt(k) + 1
    count = 0
    for n in cartesian(range(-bound, bound + 1), repeat=4):
        if sum(ni^2 for ni in n) == k:
            count += 1
    return count


# ---------------------------------------------------------------------------
# Cumulative lattice point counts
# ---------------------------------------------------------------------------

def N_d(d, K):
    """Cumulative lattice-point count: N_d(K) = |{n in Z^d : |n|^2 <= K}|."""
    bound = isqrt(K) + 1
    count = 0
    for n in cartesian(range(-bound, bound + 1), repeat=d):
        if sum(ni^2 for ni in n) <= K:
            count += 1
    return count


def N4(K):
    """N_4(K) via Jacobi: sum_{k=0}^{K} r_4(k)."""
    return sum(r4_jacobi(k) for k in range(K + 1))


# ---------------------------------------------------------------------------
# Registered tests
# ---------------------------------------------------------------------------

@register_test("number_theory", "jacobi_four_square")
def test_jacobi_four_square():
    """Verify r_4(k) = 8*sigma_tilde(k) against direct count for k=0..10."""
    assert r4_jacobi(0) == 1 == r4_direct(0)
    for k in range(1, 11):
        jacobi = r4_jacobi(k)
        direct = r4_direct(k)
        assert jacobi == direct, f"r_4({k}): Jacobi={jacobi}, direct={direct}"
    certify("jacobi_four_square", {"verified_range": [0, 10]})
    return "r_4(k) = 8*sigma_tilde(k) verified for k=0..10"


@register_test("number_theory", "cumulative_N4")
def test_cumulative_N4():
    """Verify N_4(K) values at key cutoffs."""
    expected = {1: 9, 2: 33, 3: 65, 4: 89, 5: 137}
    for K, val in expected.items():
        computed = N4(K)
        assert computed == val, f"N_4({K}) = {computed}, expected {val}"
    certify("N4_values", expected)
    return f"N_4(K) = {expected}"


@register_test("number_theory", "shell_r4")
def test_shell_r4():
    """Verify shell multiplicities r_4(k) for k=0..5."""
    expected = [1, 8, 24, 32, 24, 48]
    for k in range(6):
        assert r4_jacobi(k) == expected[k], f"r_4({k}) = {r4_jacobi(k)}, expected {expected[k]}"
    certify("shell_r4", expected)
    return f"r_4(k) = {expected} for k=0..5"


@register_test("number_theory", "dimensional_rigidity")
def test_dimensional_rigidity():
    """Verify 2d = 2^{d-1} has unique solution d=4."""
    solutions = [d for d in range(1, 100) if 2*d == 2^(d-1)]
    assert solutions == [4], f"Solutions: {solutions}"
    certify("dimensional_rigidity", {"equation": "2d = 2^{d-1}", "unique_solution": 4})
    return "2d = 2^{d-1} uniquely solved by d=4"


@register_test("number_theory", "N4_cross_verify")
def test_N4_cross_verify():
    """Cross-verify N4(K) via Jacobi formula against direct lattice count."""
    for K in range(1, 6):
        jacobi = N4(K)
        direct = N_d(4, K)
        assert jacobi == direct, f"N_4({K}): Jacobi={jacobi}, direct={direct}"
    certify("N4_cross_verified", True)
    return "N_4(K) Jacobi = N_d(4,K) direct for K=1..5"


@register_test("number_theory", "lattice_count_cross_domain")
def test_lattice_count_cross_domain():
    """Cross-verify: sum(parity_weight_counts(n,K)) from combinatorics must equal
    N_d(n,K) from number_theory. Two independent lattice-counting implementations."""
    try:
        parity_weight_counts  # from combinatorics domain
    except NameError:
        return "SKIP: combinatorics domain not loaded (standalone mode)"
    for n in [2, 3, 4]:
        for K in [n, n + 1, n + 2]:
            from_nt = N_d(n, K)
            from_comb = sum(parity_weight_counts(n, K, q=2))
            assert from_nt == from_comb, \
                f"n={n}, K={K}: N_d={from_nt}, sum(c_w)={from_comb}"
    certify("lattice_count_cross_domain", {"n_range": [2, 4], "verified": True})
    return "N_d(n,K) == sum(c_w(n,K)) verified for n=2,3,4"
