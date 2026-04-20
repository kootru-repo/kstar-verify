"""
Domain: Linear Algebra
======================
Gram matrices, spectral decomposition, rank, kernel certificates.
All computations over QQ or exact algebraic numbers. No floating point.
"""

# Framework (register_test, certify, QQ, ZZ, etc.) loaded by run_all.sage

from itertools import product as cartesian

# ---------------------------------------------------------------------------
# Gram matrix construction (lattice → association scheme)
# ---------------------------------------------------------------------------

def fixed_points_Z2(d):
    """Fixed points of Z_2 on T^d: all v in {0,1}^d."""
    return list(cartesian(range(2), repeat=d))


def hamming_dist(v1, v2):
    """Hamming distance between two binary vectors."""
    return sum((a - b) % 2 for a, b in zip(v1, v2))


def lattice_modes(d, K):
    """All n in Z^d with |n|^2 <= K."""
    bound = isqrt(K) + 1
    modes = []
    for n in cartesian(range(-bound, bound + 1), repeat=d):
        if sum(x^2 for x in n) <= K:
            modes.append(n)
    return modes


def gram_matrix(d, K):
    """Gram matrix G_{ij} = sum_n (-1)^{n.(v_i-v_j mod 2)} over |n|^2 <= K.
    Returns matrix over ZZ (exact integer)."""
    modes = lattice_modes(d, K)
    fps = fixed_points_Z2(d)
    n_fix = len(fps)
    G = matrix(ZZ, n_fix, n_fix)
    for n in modes:
        for i in range(n_fix):
            for j in range(n_fix):
                diff = tuple((fps[i][l] - fps[j][l]) % 2 for l in range(d))
                dot_mod2 = sum(n[l] * diff[l] for l in range(d)) % 2
                G[i, j] += 1 if dot_mod2 == 0 else -1
    return G


def verify_hamming_constant(G, d):
    """Verify G depends only on Hamming distance (association scheme property)."""
    fps = fixed_points_Z2(d)
    n_fix = len(fps)
    g_by_dist = {}
    for i in range(n_fix):
        for j in range(n_fix):
            h = hamming_dist(fps[i], fps[j])
            val = G[i, j]
            if h in g_by_dist:
                assert g_by_dist[h] == val, \
                    f"G not constant on Hamming dist {h}: {val} vs {g_by_dist[h]}"
            else:
                g_by_dist[h] = val
    return g_by_dist


def gram_distance_values(d, K):
    """Return G^(h) for h=0..d."""
    G = gram_matrix(d, K)
    return verify_hamming_constant(G, d)


# ---------------------------------------------------------------------------
# Spectral decomposition via Krawtchouk
# ---------------------------------------------------------------------------

def spectral_decomposition(g_vals, n, q=2):
    """Compute eigenvalues lambda_w from distance values g^(h).
    lambda_w = (1/A_w) * sum_h g(h) * K_w(h) * A_h
    where A_w = C(n,w)*(q-1)^w."""
    eigs = []
    for w in range(n + 1):
        A_w = binomial(n, w) * (q - 1)^w
        raw = QQ(0)
        for h in range(n + 1):
            A_h = binomial(n, h) * (q - 1)^h
            # Independent Krawtchouk computation (cross-checks combinatorics domain)
            K_wh = QQ(0)
            for j in range(w + 1):
                if j > h or w - j > n - h:
                    continue
                K_wh += (-1)^j * (q-1)^(w-j) * binomial(h, j) * binomial(n-h, w-j)
            raw += g_vals[h] * K_wh * A_h
        eigs.append(raw / A_w)
    return eigs


# ---------------------------------------------------------------------------
# Registered tests
# ---------------------------------------------------------------------------

@register_test("linear_algebra", "hamming_constancy")
def test_hamming_constancy():
    """Verify Gram matrix depends only on Hamming distance (association scheme).
    For each n and K, every pair (i,j) with the same Hamming distance must have
    the same G_{ij} value. This is the structural property enabling Krawtchouk
    diagonalization."""
    for n in [2, 3, 4]:
        for K in [n, n + 1]:
            G = gram_matrix(n, K)
            g_vals = verify_hamming_constant(G, n)
            # Also verify we got exactly n+1 distinct distance values
            assert len(g_vals) == n + 1, \
                f"Expected {n+1} distance values at n={n}, K={K}, got {len(g_vals)}"
    certify("hamming_constancy", {"n_range": [2, 4], "K_tested": "n and n+1"})
    return "Gram matrix Hamming-constant (association scheme) for n=2,3,4"


@register_test("linear_algebra", "distance_values_n4_K5")
def test_distance_values_n4_K5():
    """Verify exact G^(h) distance values at n=4, K=5 match paper claims.
    G^(h) = {137, 5, 17, -19, -39} for h=0..4."""
    n, K = 4, 5
    g_vals = gram_distance_values(n, K)
    expected = {0: 137, 1: 5, 2: 17, 3: -19, 4: -39}
    for h in range(n + 1):
        assert g_vals[h] == expected[h], \
            f"G^({h}) = {g_vals[h]}, expected {expected[h]}"
    certify("distance_values_n4_K5", expected)
    return f"G^(h) at n=4, K=5: {[expected[h] for h in range(n+1)]}"


@register_test("linear_algebra", "spectral_decomp_qary")
def test_spectral_decomp_qary():
    """Verify G(K) eigenvalues via Krawtchouk transform match actual matrix eigenvalues.
    NOTE: gram_matrix() constructs the Z_2 orbifold Gram, so q=2 only.
    For q>2 the Gram construction differs (Z_q orbifold) — not implemented here."""
    q = 2
    for n in [2, 3, 4]:
        K_max = 6  # enough to reach K*
        G = gram_matrix(n, K_max)
        # Get actual eigenvalues from the matrix
        actual_eigs = sorted(G.eigenvalues())
        # Get Krawtchouk-predicted eigenvalues
        g_vals_dict = verify_hamming_constant(G, n)
        g_vals = [g_vals_dict[h] for h in range(n + 1)]
        kraw_eigs = spectral_decomposition(g_vals, n, q)
        # Krawtchouk eigenvalues should match actual (with multiplicities C(n,w))
        predicted_eigs = sorted(sum(([e] * binomial(n, w)
                                     for w, e in enumerate(kraw_eigs)), []))
        assert actual_eigs == predicted_eigs, \
            f"Eigenvalue mismatch at n={n}: actual={actual_eigs}, predicted={predicted_eigs}"
    certify("spectral_decomp_verified", {"q": 2, "n_range": [2, 4]})
    return "Spectral decomposition: Krawtchouk eigenvalues match matrix eigenvalues for n=2,3,4"


def _parity_weight_counts_local(n, K, q=2):
    """Parity-weight counts c_w(K): number of lattice points m with |m|^2 <= K
    and parity weight w (number of nonzero components mod q).
    Local implementation so --domain linear_algebra works standalone."""
    bound = isqrt(K) + 1
    c_w = [0] * (n + 1)
    for m in cartesian(range(-bound, bound + 1), repeat=n):
        if sum(x^2 for x in m) <= K:
            w = sum(1 for x in m if x % q != 0)
            c_w[w] += 1
    return c_w


@register_test("linear_algebra", "parity_weight_cross_check")
def test_parity_weight_cross_check():
    """Cross-verify _parity_weight_counts_local against combinatorics domain's
    parity_weight_counts when both domains are loaded (full suite).
    Ensures the duplicated implementation is correct."""
    try:
        parity_weight_counts  # from combinatorics domain
    except NameError:
        # Running --domain linear_algebra only; skip cross-check
        return "SKIP: combinatorics domain not loaded (standalone mode)"
    for n in [2, 3, 4]:
        for K in [n, n + 1]:
            local = _parity_weight_counts_local(n, K, q=2)
            canonical = parity_weight_counts(n, K, q=2)
            assert local == canonical, \
                f"Mismatch at n={n}, K={K}: local={local}, canonical={canonical}"
    certify("parity_weight_cross_check", {"n_range": [2, 4], "verified": True})
    return "Local parity_weight_counts matches combinatorics domain for n=2,3,4"


@register_test("linear_algebra", "completeness_kn")
def test_completeness_kn():
    """Verify at K=n, all parity-weight classes are represented (c_w > 0),
    so the Gram matrix has all eigenvalues positive and is full rank."""
    for n in [2, 3, 4]:
        c_w = _parity_weight_counts_local(n, n, q=2)
        assert all(c > 0 for c in c_w), \
            f"Incomplete weight coverage at K=n={n}: c_w={c_w}"

        # Verify eigenvalues match
        g_vals_dict = gram_distance_values(n, n)
        g_vals = [g_vals_dict[h] for h in range(n + 1)]
        eigs = spectral_decomposition(g_vals, n, q=2)
        assert all(e > 0 for e in eigs), \
            f"Not all eigenvalues positive at K=n={n}: {eigs}"

        certify(f"completeness_Kn_n{n}", {"c_w": c_w, "eigs": [str(e) for e in eigs]})
    return "K=n gives all c_w>0 and all lambda_w>0 for n=2,3,4"


@register_test("linear_algebra", "rank_full")
def test_rank_full():
    """Verify rank G(K*) = 2^n for n=2..5, and rank G(K*-1) < 2^n (minimality)."""
    for n in range(2, 6):
        K_star = None
        for K in range(1, 20):
            G = gram_matrix(n, K)
            if G.rank() == 2^n:
                K_star = K
                break
        assert K_star is not None, f"Full rank not achieved for n={n} up to K=19"
        # Verify minimality: K*-1 must NOT have full rank
        if K_star >= 2:
            G_prev = gram_matrix(n, K_star - 1)
            assert G_prev.rank() < 2^n, \
                f"rank G({K_star-1}) = {G_prev.rank()} is already full rank at n={n}"
        certify(f"rank_full_n{n}", {"K_star": K_star, "rank": 2^n})
    return "rank G(K*) = 2^n and rank G(K*-1) < 2^n verified for n=2..5"


@register_test("linear_algebra", "delsarte_certificate")
def test_delsarte_certificate():
    """Verify Delsarte dual certificate: at K=3, n=4, c_4=0 so lambda_4=0,
    meaning weight-4 idempotent E_4 lies in ker G(3). At K=4, all c_w > 0
    and Gram becomes full rank."""
    n = 4
    # K=3: c_4(3) = 0 because |m|^2 <= 3 cannot have all 4 components nonzero
    g_vals_3 = gram_distance_values(n, 3)
    g_list_3 = [g_vals_3[h] for h in range(n + 1)]
    eigs_3 = spectral_decomposition(g_list_3, n)
    assert eigs_3[4] == 0, f"lambda_4(K=3) = {eigs_3[4]}, expected 0"
    assert eigs_3[3] > 0, f"lambda_3(K=3) = {eigs_3[3]}, expected > 0"

    # Verify via Gram rank: rank should be 2^4 - C(4,4) = 15 (missing weight-4 subspace)
    G3 = gram_matrix(n, 3)
    assert G3.rank() == 2^n - binomial(n, 4), \
        f"rank G(3) = {G3.rank()}, expected {2^n - binomial(n, 4)}"

    # K=4: all eigenvalues positive, full rank
    g_vals_4 = gram_distance_values(n, 4)
    g_list_4 = [g_vals_4[h] for h in range(n + 1)]
    eigs_4 = spectral_decomposition(g_list_4, n)
    assert all(e > 0 for e in eigs_4), f"K=4 eigenvalues: {eigs_4}"
    G4 = gram_matrix(n, 4)
    assert G4.rank() == 2^n, f"rank G(4) = {G4.rank()}, expected {2^n}"

    certify("delsarte_certificate_n4", {
        "K3_zero_eig_weights": [w for w in range(n+1) if eigs_3[w] == 0],
        "K3_rank": int(G3.rank()), "K4_rank": int(G4.rank()),
    })
    return f"K=3: lambda_4=0, rank={G3.rank()}; K=4: all lambda>0, rank={G4.rank()}"
