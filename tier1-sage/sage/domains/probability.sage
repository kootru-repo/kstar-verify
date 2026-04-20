"""
Domain: Probability
===================
Hypergeometric sampling, coupon-collector bounds, product inequalities.
All computations over QQ. No floating point.
"""

# Framework (register_test, certify, QQ, ZZ, etc.) loaded by run_all.sage

# ---------------------------------------------------------------------------
# Hypergeometric exact computation
# ---------------------------------------------------------------------------

def hypergeometric_exact(M, N, A):
    """Exact probability that all A special items are drawn in M draws from N.
    P = C(N-A, M-A) / C(N, M) = prod_{j=0}^{A-1} (M-j)/(N-j).
    Returns exact rational."""
    if M < A:
        return QQ(0)
    prob = QQ(1)
    for j in range(A):
        prob *= QQ(M - j) / QQ(N - j)
    return prob


def hypergeometric_upper_bound(M, N, A):
    """Upper bound: prod (M-j)/(N-j) <= (M/N)^A.
    Valid because (M-j)/(N-j) <= M/N for j >= 0, M < N."""
    return (QQ(M) / QQ(N))^A


# ---------------------------------------------------------------------------
# Registered tests
# ---------------------------------------------------------------------------

@register_test("probability", "hypergeometric_exact")
def test_hypergeometric_exact():
    """Verify P(all 12 wt-1 ops survive) in M=137 from N=255."""
    M, N, A_1 = 137, 255, 12
    p_exact = hypergeometric_exact(M, N, A_1)
    # Must be < 5e-4 (manuscript claim)
    assert p_exact < QQ(5) / QQ(10000), f"P = {float(p_exact):.6e} >= 5e-4"
    # Must be > 4e-4 (tighter than naive (M/N)^A bound)
    assert p_exact > QQ(4) / QQ(10000), f"P = {float(p_exact):.6e} <= 4e-4"
    certify("hypergeometric_w1_n4", str(p_exact), "exact_rational")
    certify("hypergeometric_w1_n4_float", float(p_exact), "float_approx")
    return f"P(all wt-1 survive) = {float(p_exact):.4e} in (4e-4, 5e-4)"


@register_test("probability", "hypergeometric_mean")
def test_hypergeometric_mean():
    """Verify E[captured from class] = A * M / N by direct enumeration.
    For small N, enumerate all C(N,M) subsets and count."""
    # Small case: N=7 items, A=3 special, draw M=4
    N_small, M_small, A_small = 7, 4, 3
    # E[captured] = A * M / N = 3 * 4 / 7 = 12/7
    formula_mean = QQ(A_small * M_small) / QQ(N_small)
    # Direct: average over all C(7,4)=35 draws
    from itertools import combinations
    total_captured = 0
    special = set(range(A_small))  # items {0,1,2} are special
    draws = list(combinations(range(N_small), M_small))
    for draw in draws:
        total_captured += len(special.intersection(draw))
    direct_mean = QQ(total_captured) / QQ(len(draws))
    assert direct_mean == formula_mean, \
        f"Direct mean {direct_mean} != formula {formula_mean}"

    # Verify the paper's weight-class sizes and their partition of operator space
    n = 4
    A_w = [binomial(n, w) * 3^w for w in range(n + 1)]
    assert A_w == [1, 12, 54, 108, 81], f"A_w = {A_w}"
    assert sum(A_w) == 4^n, f"sum(A_w) = {sum(A_w)}, expected {4^n}"
    # Non-identity operators = N = 255 (the hypergeometric pool)
    assert sum(A_w[1:]) == 4^n - 1, \
        f"Non-identity operators: {sum(A_w[1:])}, expected {4^n - 1}"
    # Linearity of expectation: sum of per-class expected captures = M
    M, N = 137, 255
    total_expected = sum(QQ(A_w[w] * M) / QQ(N) for w in range(1, n + 1))
    assert total_expected == M, \
        f"Sum of per-class E[captured] = {total_expected}, expected {M}"
    certify("hypergeometric_mean_verified", {"small_case": "N=7,M=4,A=3", "formula": "A*M/N"})
    return f"E[captured] = A*M/N verified by enumeration (N=7); sum(A_w[1:])=255, linearity check"


@register_test("probability", "hypergeometric_product")
def test_hypergeometric_product():
    """Verify prod bound: P(all A captured) <= (M/N)^A for all weight classes."""
    M, N = 137, 255
    n = 4
    A_w = [binomial(n, w) * 3^w for w in range(1, n + 1)]  # [12, 54, 108, 81]
    results = {}
    for w_idx, A in enumerate(A_w):
        w = w_idx + 1
        exact = hypergeometric_exact(M, N, A)
        bound = hypergeometric_upper_bound(M, N, A)
        assert exact <= bound, \
            f"w={w}: exact={float(exact):.4e} > bound={float(bound):.4e}"
        # Bound must be strict (exact < bound for A >= 1)
        assert exact < bound, \
            f"w={w}: bound not strict — exact == bound"
        results[f"w{w}"] = {"exact": float(exact), "bound": float(bound), "A": int(A)}
    certify("hypergeometric_product_all_weights", results)
    return f"P(all A_w captured) <= (M/N)^A_w verified for w=1..4 (A={[int(a) for a in A_w]})"


@register_test("probability", "hypergeometric_bound_value")
def test_hypergeometric_bound_value():
    """Verify exact hypergeometric bound (M/N)^A_1 = (137/255)^12.
    Registry claim: hypergeometric_bound_w1 = 5.78e-04."""
    M, N, A_1 = 137, 255, 12
    bound = (QQ(M) / QQ(N))^A_1
    # Verify it matches the registry's approximate value
    assert bound > QQ(5) / QQ(10000), f"bound = {float(bound):.4e} <= 5e-4"
    assert bound < QQ(6) / QQ(10000), f"bound = {float(bound):.4e} >= 6e-4"
    # More precisely: 5.78e-4
    assert bound > QQ(577) / QQ(1000000) and bound < QQ(579) / QQ(1000000), \
        f"bound = {float(bound):.6e}, expected ~5.78e-4"
    certify("hypergeometric_bound_w1", float(bound), "float_approx")
    certify("hypergeometric_bound_w1_exact", str(bound), "exact_rational")
    return f"(M/N)^A_1 = (137/255)^12 = {float(bound):.4e}"


@register_test("probability", "random_vanishing")
def test_random_vanishing():
    """Verify (M/N)^A_w vanishes for large A_w when M/N < 1."""
    # For n=4: M=137, N=255, ratio = 137/255
    ratio = QQ(137) / QQ(255)
    assert ratio < 1
    # A_1 = 12: bound = ratio^12
    for A in [12, 54, 108]:
        bound = ratio^A
        assert bound < QQ(1) / QQ(100), f"(M/N)^{A} = {float(bound):.4e} not small"
    certify("random_vanishing_n4", {
        "M_over_N": str(ratio),
        "bound_A12": float(ratio^12),
        "bound_A54": float(ratio^54),
    })
    return f"M/N = {float(ratio):.4f}, bounds vanish: {float(ratio^12):.2e}, {float(ratio^54):.2e}"
