"""Krawtchouk polynomial properties (Section III.C, SM S5)."""

from sympy import Matrix, binomial

from common import krawtchouk_exact, krawtchouk_matrix


def test_krawtchouk():
    """Verify Krawtchouk polynomial properties (exact rational)."""
    passed = 0
    n = 4

    K = krawtchouk_matrix(n)

    # Orthogonality relation
    weights = Matrix(n + 1, n + 1, lambda i, j:
        binomial(n, i) if i == j else 0)
    product = K * weights * K.T
    expected = Matrix(n + 1, n + 1, lambda i, j:
        2**n * binomial(n, i) if i == j else 0)
    assert product == expected, "Krawtchouk orthogonality failed"
    passed += 1
    print("  [PASS] Orthogonality: K diag(C) K^T = 2^n diag(C)")

    for h in range(5):
        assert krawtchouk_exact(0, h, 4) == 1
    passed += 1
    print("  [PASS] K_0(h;4) = 1 for all h")

    for h in range(5):
        assert krawtchouk_exact(1, h, 4) == 4 - 2 * h
    passed += 1
    print("  [PASS] K_1(h;4) = 4 - 2h")

    k2_expected = [6, 0, -2, 0, 6]
    for h in range(5):
        assert krawtchouk_exact(2, h, 4) == k2_expected[h]
    passed += 1
    print("  [PASS] K_2(h;4) = [6, 0, -2, 0, 6]")

    # Duality: K_w(h)*C(n,h) = K_h(w)*C(n,w)
    for w in range(5):
        for h in range(5):
            lhs = krawtchouk_exact(w, h, n) * binomial(n, h)
            rhs = krawtchouk_exact(h, w, n) * binomial(n, w)
            assert lhs == rhs
    passed += 1
    print("  [PASS] Krawtchouk symmetry (all 25 pairs)")

    K_full = [[krawtchouk_exact(w, h, 4) for h in range(5)]
              for w in range(5)]
    expected_matrix = [
        [1,  1,  1,  1,  1],
        [4,  2,  0, -2, -4],
        [6,  0, -2,  0,  6],
        [4, -2,  0,  2, -4],
        [1, -1,  1, -1,  1],
    ]
    assert K_full == expected_matrix
    passed += 1
    print("  [PASS] Full 5x5 Krawtchouk matrix matches manuscript")

    print(f"  [PASS] Krawtchouk polynomials (exact): {passed} checks")
    return passed
