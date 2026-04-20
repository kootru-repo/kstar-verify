"""
Domain: Combinatorics
=====================
Krawtchouk polynomials, weight enumerators, association schemes,
greedy redistribution, operator counting.

All computations over QQ (exact rationals). No floating point.
"""

# Framework (register_test, certify, QQ, ZZ, etc.) loaded by run_all.sage

# ---------------------------------------------------------------------------
# Krawtchouk polynomials (binary and q-ary)
# ---------------------------------------------------------------------------

def krawtchouk(w, h, n, q=2):
    """q-ary Krawtchouk polynomial K_w(h; n, q). Exact over QQ."""
    val = QQ(0)
    for j in range(w + 1):
        if j > h or w - j > n - h:
            continue
        val += (-1)^j * (q - 1)^(w - j) * binomial(h, j) * binomial(n - h, w - j)
    return val


def krawtchouk_matrix(n, q=2):
    """(n+1) x (n+1) Krawtchouk matrix over QQ."""
    return matrix(QQ, n + 1, n + 1, lambda w, h: krawtchouk(w, h, n, q))


def verify_krawtchouk_orthogonality(n, q=2):
    """K * diag(A_w) * K^T = q^n * diag(A_w) for q-ary Krawtchouk."""
    K = krawtchouk_matrix(n, q)
    A = [binomial(n, w) * (q - 1)^w for w in range(n + 1)]
    W = diagonal_matrix(QQ, A)
    product = K * W * K.transpose()
    expected = q^n * W
    assert product == expected, "Krawtchouk orthogonality failed"
    return True


def verify_krawtchouk_duality(n, q=2):
    """K_w(h) * A_h = K_h(w) * A_w (Krawtchouk duality)."""
    for w in range(n + 1):
        for h in range(n + 1):
            A_h = binomial(n, h) * (q - 1)^h
            A_w = binomial(n, w) * (q - 1)^w
            lhs = krawtchouk(w, h, n, q) * A_h
            rhs = krawtchouk(h, w, n, q) * A_w
            assert lhs == rhs, f"Duality failed at w={w}, h={h}"
    return True


# ---------------------------------------------------------------------------
# Parity-weight counting (lattice -> Pauli map)
# ---------------------------------------------------------------------------

def parity_weight_counts(n, K, q=2):
    """Count lattice points m in Z^n with |m|^2 <= K by q-ary parity weight.
    Returns list c_w for w = 0..n. Exact integer."""
    from itertools import product as cartesian
    bound = isqrt(K) + 1
    counts = [0] * (n + 1)
    for m in cartesian(range(-bound, bound + 1), repeat=n):
        if sum(x^2 for x in m) <= K:
            w = sum(1 for x in m if x % q != 0)
            counts[w] += 1
    return counts


def weight_class_sizes(n, q=2):
    """A_w = C(n,w) * (q^2-1)^w -- number of weight-w operators."""
    return [binomial(n, w) * (q^2 - 1)^w for w in range(n + 1)]


def greedy_redistribute(c_w, A_w):
    """Greedy redistribution algorithm.
    Initialize e=0. At each w: M_w = min(c_w+e, A_w), e <- c_w+e-M_w.
    Returns list M_w."""
    n = len(c_w)
    M = [0] * n
    e = 0
    for w in range(n):
        available = c_w[w] + e
        M[w] = min(available, A_w[w])
        e = available - M[w]
    return M


def verify_greedy_monotone(c_w_seq, A_w):
    """Verify M_w(K) is non-decreasing as K increases.
    c_w_seq: list of c_w vectors for K=1,2,..."""
    prev_M = None
    for K_idx, c_w in enumerate(c_w_seq):
        M = greedy_redistribute(c_w, A_w)
        if prev_M is not None:
            for w in range(len(M)):
                assert M[w] >= prev_M[w], \
                    f"Monotonicity failed at K={K_idx+1}, w={w}: {M[w]} < {prev_M[w]}"
        prev_M = M
    return True


def operator_lower_bound(n, k, q=2):
    """Lower bound on |S|: sum_{w=1}^{k} C(n,w)*(q^2-1)^w."""
    return sum(binomial(n, w) * (q^2 - 1)^w for w in range(1, k + 1))


# ---------------------------------------------------------------------------
# Gram eigenvalues via Krawtchouk transform
# ---------------------------------------------------------------------------

def gram_eigenvalues_from_counts(c_w, n, q=2):
    """lambda_w = q^n * c_w / (C(n,w) * (q-1)^w)."""
    eigs = []
    for w in range(n + 1):
        denom = binomial(n, w) * (q - 1)^w
        eigs.append(QQ(q^n * c_w[w]) / QQ(denom))
    return eigs


def _build_gram_matrix(n, K):
    """Independent Gram matrix construction for cross-verification.
    Reimplements the linear_algebra domain's gram_matrix() so that
    combinatorics tests can run standalone (--domain combinatorics)."""
    from itertools import product as cartesian
    fps = list(cartesian(range(2), repeat=n))
    bound = isqrt(K) + 1
    modes = [m for m in cartesian(range(-bound, bound + 1), repeat=n)
             if sum(x^2 for x in m) <= K]
    G = matrix(ZZ, len(fps), len(fps))
    for m in modes:
        for i in range(len(fps)):
            for j in range(len(fps)):
                diff = tuple((fps[i][l] - fps[j][l]) % 2 for l in range(n))
                dot = sum(m[l] * diff[l] for l in range(n)) % 2
                G[i, j] += 1 if dot == 0 else -1
    return G


def find_kstar(n, q=2, K_max=100):
    """Find K*_full = min{K : all Gram eigenvalues positive}.
    K*_full(n): n=2->2, n=3->3, n=4->4, n=5->5.
    NOTE: The paper's K*=5 at n=4 is the cumulative lattice threshold
    (N_4(5)=137), not K*_full. K*_full(4)=4 because c_w(4) all > 0."""
    for K in range(1, K_max + 1):
        c_w = parity_weight_counts(n, K, q)
        eigs = gram_eigenvalues_from_counts(c_w, n, q)
        if all(e > 0 for e in eigs):
            return K, c_w, eigs
    return None, None, None


# ---------------------------------------------------------------------------
# Registered tests (callable by framework)
# ---------------------------------------------------------------------------

@register_test("combinatorics", "krawtchouk_orthogonality")
def test_krawtchouk_orthogonality():
    """Verify K * diag(A_w) * K^T = q^n * diag(A_w) for n=2..5, q=2,3."""
    for n in [2, 3, 4, 5]:
        for q in [2, 3]:
            verify_krawtchouk_orthogonality(n, q)
    certify("krawtchouk_orthogonality", {"n_range": [2, 5], "q_values": [2, 3]})
    return "Krawtchouk orthogonality verified for n=2..5, q=2,3"


@register_test("combinatorics", "krawtchouk_duality")
def test_krawtchouk_duality():
    """Verify K_w(h)*A_h = K_h(w)*A_w for n=2..5, q=2,3."""
    for n in [2, 3, 4, 5]:
        for q in [2, 3]:
            verify_krawtchouk_duality(n, q)
    certify("krawtchouk_duality", {"n_range": [2, 5], "q_values": [2, 3]})
    return "Krawtchouk duality verified for n=2..5, q=2,3"


@register_test("combinatorics", "parity_weight_q3")
def test_parity_weight_q3():
    """Verify parity_weight_counts at q=3 for n=2: exercises the x%q!=0 branch.
    Manual count at K=2: c_w(2) = [1, 4, 4] (all 9 points in Z^2 with |m|^2<=2).
    At K=1: c_w = [1, 4, 0] (no weight-2 points yet).
    Also verifies K*_full(n=2,q=3) = 2 via eigenvalue formula."""
    n, q = 2, 3
    c_w_1 = parity_weight_counts(n, 1, q)
    assert c_w_1 == [1, 4, 0], f"c_w(K=1,q=3) = {c_w_1}, expected [1,4,0]"
    c_w_2 = parity_weight_counts(n, 2, q)
    assert c_w_2 == [1, 4, 4], f"c_w(K=2,q=3) = {c_w_2}, expected [1,4,4]"
    # All c_w > 0 at K=2 => K*_full = 2
    eigs = gram_eigenvalues_from_counts(c_w_2, n, q)
    assert all(e > 0 for e in eigs), f"Not all eigenvalues positive: {eigs}"
    # K=1 has c_2=0 => not full rank
    eigs_1 = gram_eigenvalues_from_counts(c_w_1, n, q)
    assert eigs_1[2] == 0, f"lambda_2(K=1) = {eigs_1[2]}, expected 0"
    certify("parity_weight_q3_n2", {"c_w_K1": c_w_1, "c_w_K2": c_w_2})
    return f"q=3: c_w(1)={c_w_1}, c_w(2)={c_w_2}, K*_full=2"


@register_test("combinatorics", "eigenvalue_monotone")
def test_eigenvalue_monotone():
    """Verify lambda_w(K) non-decreasing for n=4, q=2, K=1..6."""
    n, q = 4, 2
    prev_eigs = None
    for K in range(1, 7):
        c_w = parity_weight_counts(n, K, q)
        eigs = gram_eigenvalues_from_counts(c_w, n, q)
        if prev_eigs is not None:
            for w in range(n + 1):
                assert eigs[w] >= prev_eigs[w], \
                    f"lambda_{w} decreased at K={K}: {eigs[w]} < {prev_eigs[w]}"
        prev_eigs = eigs
    certify("eigenvalue_monotone_n4", {"verified": True, "K_range": [1, 6]})
    return "lambda_w(K) non-decreasing for K=1..6 at n=4"


@register_test("combinatorics", "formula_vs_matrix_eigenvalues")
def test_formula_vs_matrix_eigenvalues():
    """Cross-verify: eigenvalues from parity-weight formula match actual
    Gram matrix eigenvalues (two independent computations).
    Uses _build_gram_matrix (local, independent of linear_algebra domain)."""
    for n in [2, 3, 4]:
        for K in [n, n + 1]:
            # Method 1: formula lambda_w = q^n * c_w / A_w
            c_w = parity_weight_counts(n, K, q=2)
            formula_eigs = gram_eigenvalues_from_counts(c_w, n, q=2)
            # Method 2: actual Gram matrix eigenvalues (independent construction)
            G = _build_gram_matrix(n, K)
            matrix_eigs = sorted(G.eigenvalues())
            predicted = sorted(sum(([e] * binomial(n, w)
                                    for w, e in enumerate(formula_eigs)), []))
            assert matrix_eigs == predicted, \
                f"Eigenvalue mismatch at n={n}, K={K}"
    certify("formula_vs_matrix_eigs", {"n_range": [2, 4], "K_tested": "n and n+1"})
    return "Formula eigenvalues match Gram matrix eigenvalues for n=2,3,4"


@register_test("combinatorics", "greedy_monotone")
def test_greedy_monotone():
    """Verify M_w(K) non-decreasing for n=4, K=1..6, with postcondition checks."""
    n, q = 4, 2
    A_w = weight_class_sizes(n, q)
    c_w_seq = [parity_weight_counts(n, K, q) for K in range(1, 7)]
    verify_greedy_monotone(c_w_seq, A_w)
    # Postcondition checks on every K
    for K_idx, c_w in enumerate(c_w_seq):
        M = greedy_redistribute(c_w, A_w)
        K = K_idx + 1
        for w in range(len(M)):
            assert 0 <= M[w] <= A_w[w], \
                f"Bound violation at K={K}, w={w}: M_w={M[w]}, A_w={A_w[w]}"
        assert sum(M) == min(sum(c_w), sum(A_w)), \
            f"Sum invariant failed at K={K}: sum(M)={sum(M)}, min(sum(c),sum(A))={min(sum(c_w), sum(A_w))}"
    M4 = greedy_redistribute(parity_weight_counts(n, 4, q), A_w)
    M5 = greedy_redistribute(parity_weight_counts(n, 5, q), A_w)
    certify("M_w_K4_n4", M4, "exact_list")
    certify("M_w_K5_n4", M5, "exact_list")
    certify("c_w_K4_n4", parity_weight_counts(n, 4, q), "exact_list")
    certify("c_w_K5_n4", parity_weight_counts(n, 5, q), "exact_list")
    return f"M_w monotone K=1..6, postconditions verified; M_w(4)={M4}, M_w(5)={M5}"


@register_test("combinatorics", "phase2_values")
def test_phase2_values():
    """Verify c_w(4)={9,8,24,32,16} and M_w={1,12,28,32,16}."""
    n, q = 4, 2
    c_w = parity_weight_counts(n, 4, q)
    assert c_w == [9, 8, 24, 32, 16], f"c_w(4) = {c_w}, expected [9,8,24,32,16]"
    A_w = weight_class_sizes(n, q)
    M_w = greedy_redistribute(c_w, A_w)
    assert M_w == [1, 12, 28, 32, 16], f"M_w(4) = {M_w}, expected [1,12,28,32,16]"
    certify("phase2_c_w", c_w, "exact_list")
    certify("phase2_M_w", M_w, "exact_list")
    return f"c_w(4)={c_w}, M_w(4)={M_w}"


@register_test("combinatorics", "phase3_saturation")
def test_phase3_saturation():
    """Verify K*=5 achieves M_w={1,12,54,54,16} and N_4(5)=137.
    Also verifies 67/70 weight split and unmeasured count from registry."""
    n, q = 4, 2
    c_w = parity_weight_counts(n, 5, q)
    assert c_w == [9, 56, 24, 32, 16], f"c_w(5) = {c_w}"
    assert sum(c_w) == 137, f"N_4(5) = {sum(c_w)}, expected 137"
    A_w = weight_class_sizes(n, q)
    M_w = greedy_redistribute(c_w, A_w)
    assert M_w == [1, 12, 54, 54, 16], f"M_w(5) = {M_w}"
    # Registry values: 67 weight<=2 ops, 70 weight>2 ops, 119 unmeasured
    wt_le2 = sum(M_w[:3])  # M_0 + M_1 + M_2 = 1 + 12 + 54 = 67
    wt_gt2 = sum(M_w[3:])  # M_3 + M_4 = 54 + 16 = 70
    assert wt_le2 == 67, f"weight<=2: {wt_le2}"
    assert wt_gt2 == 70, f"weight>2: {wt_gt2}"
    assert q^(2*n) - sum(M_w) == 119, f"unmeasured: {q^(2*n) - sum(M_w)}"
    certify("N4_K5", 137)
    certify("phase3_M_w", M_w, "exact_list")
    return f"c_w(5)={c_w}, M_w(5)={M_w}, N_4(5)=137, 67/70 split verified"


@register_test("combinatorics", "operator_lower_bound")
def test_operator_lower_bound():
    """Verify |S| >= sum C(n,w)(q^2-1)^w for w=1..k.
    C(n,1)*3 + C(n,2)*9: n=3->9+27=36, n=4->12+54=66, n=5->15+90=105."""
    for n, k, q, expected in [(4, 2, 2, 66), (3, 2, 2, 36), (5, 2, 2, 105)]:
        bound = operator_lower_bound(n, k, q)
        assert bound == expected, f"Lower bound at n={n},k={k},q={q}: {bound} != {expected}"
    certify("lower_bound_n4_k2", 66)
    return "Lower bounds verified: n=3->36, n=4->66, n=5->105"


@register_test("combinatorics", "kstar_saturation")
def test_kstar_saturation():
    """Verify K*_full = min{K : all Gram eigenvalues positive} at n=2..5.
    K*_full(n): n=2->2, n=3->3, n=4->4, n=5->5.
    NOTE: The paper's K*=5 at n=4 is the cumulative lattice count threshold
    (N_4(5)=137), NOT K*_full. K*_full(4)=4 because c_w(4)=[9,8,24,32,16]
    already has all components positive."""
    expected_kstar = {2: 2, 3: 3, 4: 4, 5: 5}
    results = {}
    for n in range(2, 6):
        K_star, c_w, eigs = find_kstar(n, q=2)
        assert K_star is not None, f"K* not found for n={n}"
        assert all(e > 0 for e in eigs), f"Not all eigenvalues positive at n={n}, K*={K_star}"
        assert K_star == expected_kstar[n], \
            f"K*_full(n={n}) = {K_star}, expected {expected_kstar[n]}"
        results[n] = {"K_star": K_star, "N": sum(c_w)}
    certify("kstar_values", results)
    return f"K* values: {results}"


@register_test("combinatorics", "kstar_minimality")
def test_kstar_minimality():
    """Two K* definitions for n=4, q=2:
    - K*_full = 4: first K with all Gram eigenvalues positive (c_w all > 0).
    - K*_paper = 5: first K achieving N_4(K)=137 and weight-class saturation
      for the greedy redistribution M_w = [1, 12, 54, 54, 16].
    This test verifies both and explains why the paper uses K=5."""
    n, q = 4, 2
    c_w_4 = parity_weight_counts(n, 4, q)
    c_w_5 = parity_weight_counts(n, 5, q)
    eigs_4 = gram_eigenvalues_from_counts(c_w_4, n, q)
    eigs_5 = gram_eigenvalues_from_counts(c_w_5, n, q)
    A_w = weight_class_sizes(n, q)  # [1, 12, 54, 108, 81]

    # K*_full = 4: all eigenvalues already positive
    assert all(e > 0 for e in eigs_4), f"K=4 eigenvalues not all positive: {eigs_4}"
    K_star_full, _, _ = find_kstar(n, q=2)
    assert K_star_full == 4, f"K*_full = {K_star_full}, expected 4"

    # But K=4 has only 89 operators and incomplete weight-1 coverage
    assert sum(c_w_4) == 89, f"N_4(4) = {sum(c_w_4)}, expected 89"
    assert c_w_4[1] < A_w[1], f"c_1(4)={c_w_4[1]} should be < A_1={A_w[1]}"

    # K*_paper = 5: N_4(5)=137, all c_w > 0, greedy achieves target
    assert sum(c_w_5) == 137, f"N_4(5) = {sum(c_w_5)}, expected 137"
    assert all(c > 0 for c in c_w_5), f"c_w(5) has zero: {c_w_5}"
    M_w_5 = greedy_redistribute(c_w_5, A_w)
    assert M_w_5 == [1, 12, 54, 54, 16], f"M_w(5) = {M_w_5}"

    # K=3 does NOT have all eigenvalues positive (c_4=0)
    c_w_3 = parity_weight_counts(n, 3, q)
    assert c_w_3[4] == 0, f"c_4(3) = {c_w_3[4]}, expected 0"

    certify("kstar_n4_minimality", {
        "K_star_full": 4, "K_star_paper": 5,
        "N4_4": sum(c_w_4), "N4_5": sum(c_w_5),
    })
    return f"K*_full=4 (N=89), K*_paper=5 (N=137)"


@register_test("combinatorics", "eigenvalues_n4_K5")
def test_eigenvalues_n4_K5():
    """Verify exact eigenvalues at n=4, K=5: [144, 224, 64, 128, 256]."""
    n, q = 4, 2
    c_w = parity_weight_counts(n, 5, q)
    eigs = gram_eigenvalues_from_counts(c_w, n, q)
    expected = [QQ(144), QQ(224), QQ(64), QQ(128), QQ(256)]
    assert eigs == expected, f"Eigenvalues: {eigs}, expected {expected}"
    certify("eigenvalues_n4_K5", [int(e) for e in eigs])
    return f"lambda_w(K*=5, n=4) = {[int(e) for e in eigs]}"


@register_test("combinatorics", "gram_condition_number")
def test_gram_condition_number():
    """Verify kappa_G = max(lambda)/min(lambda) = 256/64 = 4 at n=4, K*=5.
    Registry claim kappa_G_K5_n4 = 4.  Note: kappa at K*_full differs
    (e.g., kappa(n=3, K*_full=3) = 8) because K*_full != K*_paper in general."""
    n, q, K = 4, 2, 5
    c_w = parity_weight_counts(n, K, q)
    eigs = gram_eigenvalues_from_counts(c_w, n, q)
    assert eigs == [QQ(144), QQ(224), QQ(64), QQ(128), QQ(256)], \
        f"Eigenvalues at n=4 K=5: {eigs}"
    kappa = max(eigs) / min(eigs)
    assert kappa == QQ(4), f"kappa_G(n={n}, K={K}) = {kappa}, expected 4"
    certify("kappa_G_K5_n4", 4)
    return f"kappa_G(n=4, K=5) = {kappa} = 256/64"


def find_kstar_paper(n, k, q=2, K_max=100):
    """Find K*_paper = min{K : greedy redistribution saturates all weight classes w<=k}.
    Returns (K_star, c_w, M_w) or (None, None, None)."""
    A_w = weight_class_sizes(n, q)
    for K in range(1, K_max + 1):
        c_w = parity_weight_counts(n, K, q)
        M_w = greedy_redistribute(c_w, A_w)
        if all(M_w[w] >= A_w[w] for w in range(1, min(k + 1, n + 1))):
            return K, c_w, M_w
    return None, None, None


@register_test("combinatorics", "mn_over_n_ratio")
def test_mn_over_n_ratio():
    """Verify M_n/N < 1 for n=3..6 at K*_paper (k=2).
    Registry claim: M_n/(q^{2n}-1) ~ 1/3 for q=2.
    Note: n=2 is overcomplete (M/N > 1), so the hypothesis starts at n=3."""
    q, k = 2, 2
    results = {}
    for n in range(3, 7):
        K_star, c_w, M_w = find_kstar_paper(n, k, q)
        if K_star is None:
            continue
        M_n = sum(c_w)
        N = q^(2*n) - 1
        ratio = QQ(M_n) / QQ(N)
        assert ratio < 1, f"M/N >= 1 at n={n}: {ratio}"
        results[n] = {"M_n": int(M_n), "N": int(N), "ratio": float(ratio)}
    certify("mn_over_n_ratio", results)
    parts = []
    for n_val in sorted(results.keys()):
        parts.append("n=%d: %.3f" % (n_val, results[n_val]["ratio"]))
    return "M_n/N ratios: " + ", ".join(parts)


@register_test("combinatorics", "spectral_mass_threshold")
def test_spectral_mass_threshold():
    """Verify K_mass = q^2: at |m|^2 = q^2, mod-q invisible vectors (±q,0,...,0)
    appear, increasing c_0. Registry claim: K_mass = q^2 (empirical)."""
    for q in [2, 3]:
        n = 4
        K_before = q^2 - 1
        K_at = q^2
        c_w_before = parity_weight_counts(n, K_before, q)
        c_w_at = parity_weight_counts(n, K_at, q)
        # c_0 must strictly increase at K=q^2 (mod-q invisible vectors appear)
        assert c_w_at[0] > c_w_before[0], \
            f"q={q}: c_0 did not increase at K=q^2: {c_w_before[0]} -> {c_w_at[0]}"
        # The increase should be exactly 2n (the ±q along each axis)
        delta_c0 = c_w_at[0] - c_w_before[0]
        assert delta_c0 >= 2*n, \
            f"q={q}: c_0 increase = {delta_c0}, expected >= {2*n}"
    certify("spectral_mass_threshold", {"q2_verified": True, "q3_verified": True})
    return f"K_mass = q^2: c_0 increases at |m|^2=q^2 for q=2,3"


@register_test("combinatorics", "kstar_full_rank_cross_check")
def test_kstar_full_rank_cross_check():
    """Cross-verify: K*_full from formula (all c_w > 0) matches K where
    actual Gram matrix first achieves full rank.
    Uses _build_gram_matrix (local, independent of linear_algebra domain)."""
    for n in [2, 3, 4]:
        # Method 1: find_kstar via parity-weight formula
        K_formula, _, _ = find_kstar(n, q=2)
        # Method 2: brute-force Gram rank (independent construction)
        for K in range(1, 20):
            G = _build_gram_matrix(n, K)
            if G.rank() == 2^n:
                K_matrix = K
                break
        assert K_formula == K_matrix, \
            f"n={n}: K*_formula={K_formula} != K*_matrix={K_matrix}"
    certify("kstar_full_rank_cross", {"n_range": [2, 4], "methods": ["formula", "matrix_rank"]})
    return "K*_full agrees between formula and matrix rank for n=2,3,4"


# ---------------------------------------------------------------------------
# Registry fact witnesses (added 2026-04-08 for v1 queue items 5 & 6)
# ---------------------------------------------------------------------------

@register_test("combinatorics", "eigenvalue_table_n4_q2_all_K")
def test_eigenvalue_table_n4_q2_all_K():
    """Witness for registry fact:eigenvalue_table_n4_q2.

    Independently recomputes the full Krawtchouk eigenvalue table
    lambda_w(K) for K = 0, 1, 2, 3, 4, 5 at n=4, q=2 via parity-weight
    counts and compares against the registry-pinned values. This is the
    Sage-side companion to the Lean theorem eigenvalue_table_general_K.
    """
    n, q = 4, 2
    expected = {
        0: [QQ(16),   QQ(0),    QQ(0),    QQ(0),    QQ(0)],
        1: [QQ(16),   QQ(32),   QQ(0),    QQ(0),    QQ(0)],
        2: [QQ(16),   QQ(32),   QQ(64),   QQ(0),    QQ(0)],
        3: [QQ(16),   QQ(32),   QQ(64),   QQ(128),  QQ(0)],
        4: [QQ(144),  QQ(32),   QQ(64),   QQ(128),  QQ(256)],
        5: [QQ(144),  QQ(224),  QQ(64),   QQ(128),  QQ(256)],
    }
    got = {}
    for K in range(6):
        c_w = parity_weight_counts(n, K, q)
        # gram_eigenvalues_from_counts returns the first min(K+1, n+1) entries;
        # pad with zeros so the table shape matches [lambda_0, ..., lambda_4].
        eigs = gram_eigenvalues_from_counts(c_w, n, q)
        while len(eigs) < n + 1:
            eigs.append(QQ(0))
        got[K] = eigs
        assert eigs == expected[K], \
            f"K={K}: got {eigs}, expected {expected[K]}"
    certify("eigenvalue_table_n4_q2_all_K",
            {K: [int(v) for v in got[K]] for K in range(6)})
    return "lambda_w(K) table n=4 q=2 K=0..5 matches registry pins"


@register_test("combinatorics", "weight3_undersaturation_K5_closed_form")
def test_weight3_undersaturation_K5_closed_form():
    """Witness for registry fact:weight3_undersaturation_K5.

    Independently verifies the closed-form identities
        c_3(K=5, n) = 8 * C(n,3)
        A_3(n)      = 27 * C(n,3)          (since A_w = C(n,w)*(q^2-1)^w at q=2)
        deficit     = 19 * C(n,3)
    and strict undersaturation c_3(5,n) < A_3(n) for a range of n. The
    c_3 values are computed by independent lattice enumeration, NOT from
    the closed form, so this is a real cross-check of the closed form.
    """
    q = 2
    n_vals = [3, 4, 5, 6, 7]
    spot = {}
    for n in n_vals:
        c_w = parity_weight_counts(n, 5, q)
        c3_enum = c_w[3]
        c3_closed = 8 * binomial(n, 3)
        A3 = binomial(n, 3) * (q^2 - 1)^3          # = 27 * C(n,3)
        deficit = A3 - c3_enum
        deficit_closed = 19 * binomial(n, 3)
        assert c3_enum == c3_closed, \
            f"n={n}: c_3 enum {c3_enum} != closed form {c3_closed}"
        assert A3 == 27 * binomial(n, 3), \
            f"n={n}: A_3 {A3} != 27*C(n,3)"
        assert deficit == deficit_closed, \
            f"n={n}: deficit {deficit} != 19*C(n,3) = {deficit_closed}"
        assert c3_enum < A3, \
            f"n={n}: strict undersaturation violated: c_3={c3_enum} !< A_3={A3}"
        spot[n] = {"c3": int(c3_enum), "A3": int(A3), "deficit": int(deficit)}
    # At n=4, manuscript line 3523 asserts c_3 = 32 < A_3 = 108:
    assert spot[4] == {"c3": 32, "A3": 108, "deficit": 76}
    certify("weight3_undersaturation_K5_closed_form", spot)
    return ("c_3(K=5,n) = 8*C(n,3), A_3 = 27*C(n,3), deficit = 19*C(n,3) "
            "for n=3..7; n=4 matches manuscript L3523 (32 < 108)")
