"""Gram-Krawtchouk diagonalisation (Lemma 6, SM Proposition 1)."""

from sympy import Rational, binomial

from common import (
    N4, gram_matrix_exact, fixed_points_Z2,
    krawtchouk_exact, hamming_dist, binary_parity_weight_counts
)
from registry import claims


def _compute_eigenvalues(g_vals, n=4):
    """lambda_w = sum_h g(h) K_w(h) C(n,h) / C(n,w)."""
    eigenvalues = []
    for w in range(n + 1):
        raw = Rational(0)
        for h in range(n + 1):
            raw += g_vals[h] * krawtchouk_exact(w, h, n) * binomial(n, h)
        q, r = divmod(raw, binomial(n, w))
        assert r == 0, f"raw_eig[{w}] not divisible by C({n},{w})"
        eigenvalues.append(q)
    return eigenvalues


def test_gram_krawtchouk_diagonalisation():
    """The Gram matrix is diagonalised by the Krawtchouk transform."""
    passed = 0
    d = 4
    fps = fixed_points_Z2(d)

    G_full = gram_matrix_exact(d, 5, dynamical=False)

    # G must be constant on each Hamming shell (association scheme property)
    g_by_dist = {}
    for i in range(16):
        for j in range(16):
            h = hamming_dist(fps[i], fps[j])
            val = G_full[i, j]
            if h in g_by_dist:
                assert g_by_dist[h] == val, \
                    f"G not constant on Hamming distance {h}"
            else:
                g_by_dist[h] = val
    passed += 1
    print("  [PASS] Gram matrix depends only on Hamming distance")

    g_vals = [g_by_dist[h] for h in range(5)]
    eigenvalues = _compute_eigenvalues(g_vals)

    for w in range(5):
        assert eigenvalues[w] > 0
    passed += 1
    print(f"  [PASS] All eigenvalues positive at K*=5: "
          f"{[int(e) for e in eigenvalues]}")

    # lambda_w = 16 c_w / C(4,w) for q=2
    c_w = binary_parity_weight_counts(4, 5)
    for w in range(5):
        expected = Rational(16 * c_w[w], binomial(4, w))
        assert eigenvalues[w] == expected, \
            f"lambda_{w} = {eigenvalues[w]}, formula gives {expected}"
    passed += 1
    print("  [PASS] Eigenvalue formula lambda_w = 16 * c_w / C(4,w)")

    n4_M = claims.get("thm:basin", "n4_M")
    trace = sum(binomial(4, w) * eigenvalues[w] for w in range(5))
    assert trace == 16 * n4_M
    passed += 1
    print(f"  [PASS] Trace = 16 * {n4_M} = {trace}")

    expected_cw = claims.get("lem:monotone", "c_w_K5")
    assert c_w == expected_cw
    assert sum(c_w) == n4_M
    assert N4(5) == n4_M
    passed += 1
    print(f"  [PASS] c_w = {c_w}, N_4(5) = {sum(c_w)} = {n4_M}")

    assert all(c > 0 for c in c_w)
    passed += 1
    print("  [PASS] All c_w > 0 (support-completeness at K*=5)")

    print(f"  [PASS] Gram-Krawtchouk diagonalisation: {passed} checks")
    return passed
