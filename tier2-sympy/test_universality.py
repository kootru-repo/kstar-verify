"""Universality across Hamming schemes H(n,q) (Lemma 5, SM Proposition 1)."""

from sympy import Rational, binomial, Matrix

from common import (
    krawtchouk_qary, krawtchouk_qary_matrix,
    qary_parity_weight_counts, qary_gram_eigenvalues,
)
from registry import claims


def test_universality():
    """Verify universality of Krawtchouk spectral correspondence."""
    passed = 0

    # K* = q^2: c_0 jumps from 1 (only origin) to > 1 (origin + nearest shell)
    # At n=4, the nearest shell has 8 vectors -> c_0 jumps to 9.
    # This is dimension-specific, so we test the GENERAL property (jump > 1)
    # and the n=4 specific value separately.
    for q in [2, 3, 5, 7]:
        n = 4
        K_star = q**2

        c_below = qary_parity_weight_counts(n, K_star - 1, q)
        assert c_below[0] == 1, \
            f"c_0(K*-1={K_star-1}, q={q}) = {c_below[0]}, expected 1"

        c_at = qary_parity_weight_counts(n, K_star, q)
        assert c_at[0] > 1, \
            f"c_0(K*={K_star}, q={q}) = {c_at[0]}, expected > 1"

        # n=4: nearest shell in Z^4 has 2*4=8 vectors, so c_0 = 1+8 = 9
        expected_c0 = 1 + 2 * n  # origin + nearest-neighbor shell
        assert c_at[0] == expected_c0, \
            f"c_0(K*={K_star}, q={q}) = {c_at[0]}, expected {expected_c0}"

        passed += 1
        print(f"  [PASS] K*=q^2={K_star} (q={q}): "
              f"c_0 jumps from {c_below[0]} to {c_at[0]}")

    # --- q-ary orthogonality ---

    for q, n in [(2, 4), (3, 2), (5, 2)]:
        K = krawtchouk_qary_matrix(n, q)
        D_h = Matrix(n + 1, n + 1, lambda i, j:
            binomial(n, i) * (q - 1)**i if i == j else 0)
        product = K * D_h * K.T
        expected = Matrix(n + 1, n + 1, lambda i, j:
            q**n * binomial(n, i) * (q - 1)**i if i == j else 0)
        assert product == expected, \
            f"Krawtchouk orthogonality failed for q={q}, n={n}"
        passed += 1
        print(f"  [PASS] q-ary Krawtchouk orthogonality (q={q}, n={n})")

    # --- positivity at operational K* ---

    op_cases = [
        (2, 4, 5),
        (3, 4, 9),
        (5, 4, 25),
        (7, 4, 49),
    ]

    for q, n, kstar in op_cases:
        eigs = qary_gram_eigenvalues(n, kstar, q)
        assert all(e > 0 for e in eigs), \
            f"Not all eigenvalues positive at K*={kstar} (q={q}): {eigs}"
        passed += 1
        print(f"  [PASS] All eigenvalues positive at K*={kstar} "
              f"(q={q}, n={n}): {[float(e) for e in eigs]}")

    # --- operator counts ---

    n4_M = claims.get("thm:basin", "n4_M")
    c_q2 = qary_parity_weight_counts(4, 5, 2)
    assert sum(c_q2) == n4_M
    passed += 1
    print(f"  [PASS] N_4(5, q=2) = {sum(c_q2)} = {n4_M}")

    c_q3_n2 = qary_parity_weight_counts(2, 9, 3)
    assert sum(c_q3_n2) == 29, f"N_2(9, q=3) = {sum(c_q3_n2)}, expected 29"
    passed += 1
    print(f"  [PASS] N_2(9, q=3) = {sum(c_q3_n2)} = 29, c_w = {c_q3_n2}")

    # trace: sum_w C(n,w)(q-1)^w lambda_w = q^n N
    for q, n, kstar in op_cases:
        eigs = qary_gram_eigenvalues(n, kstar, q)
        c = qary_parity_weight_counts(n, kstar, q)
        N_total = sum(c)

        trace = sum(binomial(n, w) * (q - 1)**w * eigs[w]
                    for w in range(n + 1))
        expected_trace = q**n * N_total
        assert trace == expected_trace, \
            f"Trace mismatch: {trace} vs {expected_trace} (q={q})"
        passed += 1
        print(f"  [PASS] Trace identity (q={q}, n={n}): "
              f"q^n * N = {expected_trace}")

    # --- q-ary symmetry ---

    for q, n in [(3, 2), (5, 2), (3, 4)]:
        for w in range(n + 1):
            for h in range(n + 1):
                lhs = (krawtchouk_qary(w, h, n, q)
                       * binomial(n, h) * (q - 1)**h)
                rhs = (krawtchouk_qary(h, w, n, q)
                       * binomial(n, w) * (q - 1)**w)
                assert lhs == rhs, \
                    f"Symmetry failed at q={q}, w={w}, h={h}"
        passed += 1
        print(f"  [PASS] q-ary Krawtchouk symmetry (q={q}, n={n})")

    # cross-check against Gram-Krawtchouk test
    eigs_q2 = qary_gram_eigenvalues(4, 5, 2)
    expected_q2 = [Rational(e) for e in claims.get("prop:spectral_q_main", "eigenvalues_n4_K5_q2")]
    for w in range(5):
        assert eigs_q2[w] == expected_q2[w]
    passed += 1
    print(f"  [PASS] q=2 eigenvalues match: "
          f"{[int(e) for e in eigs_q2]}")

    # --- MULTI-n SWEEP: eigenvalue positivity and support-completeness ---
    # Lemma 5 is stated for general n. Verify at n=2,3,4,5,6 (q=2).
    # For each n, find K* (minimal K where all eigenvalues > 0) and verify
    # support-completeness (all c_w > 0).

    from common import weight_class_sizes

    for n_test in [2, 3, 5, 6]:
        # Find K* by scanning K = 1, 2, ... until all eigenvalues positive
        K_star_found = None
        for K_try in range(1, n_test**2 + 2):
            eigs = qary_gram_eigenvalues(n_test, K_try, 2)
            if all(e > 0 for e in eigs):
                K_star_found = K_try
                break

        assert K_star_found is not None, \
            f"n={n_test}: no K* found up to K={n_test**2+1}"

        # Verify support-completeness at K*
        c_w = qary_parity_weight_counts(n_test, K_star_found, 2)
        assert all(c > 0 for c in c_w), \
            f"n={n_test}, K*={K_star_found}: c_w has zeros: {c_w}"

        N_total = sum(c_w)

        # Verify Delsarte certificate: K*-1 has a zero eigenvalue
        if K_star_found > 1:
            eigs_below = qary_gram_eigenvalues(n_test, K_star_found - 1, 2)
            has_zero = any(e == 0 for e in eigs_below)
            assert has_zero, \
                f"n={n_test}: K*-1={K_star_found-1} should have zero eigenvalue"

        # Verify trace identity
        eigs_star = qary_gram_eigenvalues(n_test, K_star_found, 2)
        trace = sum(binomial(n_test, w) * eigs_star[w]
                    for w in range(n_test + 1))
        assert trace == 2**n_test * N_total

        passed += 1
        print(f"  [PASS] n={n_test}: K*={K_star_found}, "
              f"N={N_total}, c_w={c_w}, all c_w > 0, "
              f"trace={int(trace)}")

    print(f"  [PASS] Universality: {passed} checks")
    return passed
