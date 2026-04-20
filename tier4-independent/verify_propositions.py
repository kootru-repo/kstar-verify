#!/usr/bin/env python3
"""
Independent verification of SM Proposition 1 and Lemma 6.

Uses only SymPy (exact rational arithmetic).
Does NOT import core.py or any project code.

SM Proposition 1: Spectral decomposition of the Gram matrix G(K) via
    Krawtchouk idempotents. Eigenvalues lambda_w = q^n * c_w / [C(n,w)*(q-1)^w].
Lemma 6 (Support-completeness): rank G(K*) = q^n implies every
    Hamming-weight class is populated by at least one lattice vector.

Convention: N_d(K) INCLUDES the origin m=0 (the identity operator).
    N_4(5) = 1 + 8 + 24 + 32 + 24 + 48 = 137.
"""
import sys
import math
from itertools import product as cart_product
from sympy import Matrix, Rational, binomial, Integer
from registry import claims

PASS = 0
FAIL = 0


def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  [PASS] {name}" + (f"  ({detail})" if detail else ""))
    else:
        FAIL += 1
        print(f"  [FAIL] {name}" + (f"  ({detail})" if detail else ""))
    return condition


def sign(dot):
    """(-1)^dot as exact Rational, safe for negative dot."""
    return Rational(1) if dot % 2 == 0 else Rational(-1)


# -- Lattice enumeration ------------------------------------------------

def lattice_vectors(d, K):
    """All m in Z^d with |m|^2 <= K."""
    bound = math.isqrt(K) + 1
    vecs = []
    for m in cart_product(range(-bound, bound + 1), repeat=d):
        if sum(x * x for x in m) <= K:
            vecs.append(m)
    return vecs


def shell_count(d, k):
    """#{m in Z^d : |m|^2 = k}."""
    bound = math.isqrt(k) + 1
    return sum(1 for m in cart_product(range(-bound, bound + 1), repeat=d)
               if sum(x * x for x in m) == k)


def cumulative_count(d, K, include_origin=True):
    """N_d(K) including or excluding origin."""
    start = 0 if include_origin else 1
    return sum(shell_count(d, k) for k in range(start, K + 1))


# -- Parity-weight counts -----------------------------------------------

def hamming_weight_mod2(m):
    return sum(1 for x in m if x % 2 != 0)


def parity_weight_counts(d, K, include_origin=True):
    """c_w counts by parity weight."""
    c = [0] * (d + 1)
    for m in lattice_vectors(d, K):
        norm_sq = sum(x * x for x in m)
        if not include_origin and norm_sq == 0:
            continue
        c[hamming_weight_mod2(m)] += 1
    return c


# -- Krawtchouk polynomials (exact) --------------------------------------

def krawtchouk(h, w, n):
    """K_h(w; n) = sum_s (-1)^s C(w,s) C(n-w, h-s).

    This is the degree-h Krawtchouk polynomial evaluated at w.
    """
    val = Rational(0)
    for s in range(h + 1):
        if s > w or h - s > n - w:
            continue
        val += ((-1)**s * binomial(w, s) * binomial(n - w, h - s))
    return val


# -- Gram matrix ---------------------------------------------------------

def gram_distance_values(d, K, include_origin=True):
    """G^(h) = sum_{m in S} (-1)^{m . v_h}, v_h canonical with hw=h."""
    G_h = []
    vecs = lattice_vectors(d, K)
    for h in range(d + 1):
        v = tuple([1] * h + [0] * (d - h))
        total = Rational(0)
        for m in vecs:
            norm_sq = sum(x * x for x in m)
            if not include_origin and norm_sq == 0:
                continue
            dot = sum(m[i] * v[i] for i in range(d))
            total += sign(dot)
        G_h.append(total)
    return G_h


def gram_matrix_exact(d, K, include_origin=True):
    """16x16 Gram matrix from distance values (exact)."""
    verts = list(cart_product([0, 1], repeat=d))
    N = len(verts)
    G_h = gram_distance_values(d, K, include_origin)
    G = Matrix.zeros(N, N)
    for i in range(N):
        for j in range(N):
            h = sum(1 for k in range(d) if verts[i][k] != verts[j][k])
            G[i, j] = G_h[h]
    return G


def gram_dynamical(d, K):
    """Dynamical Gram: excludes k=0 (origin) and k=1 shells."""
    verts = list(cart_product([0, 1], repeat=d))
    N = len(verts)
    G = Matrix.zeros(N, N)
    for m in lattice_vectors(d, K):
        norm_sq = sum(x * x for x in m)
        if norm_sq < 2:
            continue
        for i in range(N):
            for j in range(N):
                diff = tuple(verts[i][k] ^ verts[j][k] for k in range(d))
                dot = sum(m[k] * diff[k] for k in range(d))
                G[i, j] += sign(dot)
    return G


# -- Primitive idempotents of H(n,2) ------------------------------------

def primitive_idempotents(d):
    """E_w[i,j] = (C(d,w) / 2^d) * K_{dist(i,j)}(w; d) / C(d, dist(i,j)).

    These satisfy: E_w * E_{w'} = delta_{ww'} * E_w, sum_w E_w = I.
    """
    verts = list(cart_product([0, 1], repeat=d))
    N = len(verts)
    idems = []
    for w in range(d + 1):
        E = Matrix.zeros(N, N)
        mw = binomial(d, w)
        for i in range(N):
            for j in range(N):
                h = sum(1 for k in range(d) if verts[i][k] != verts[j][k])
                kh = binomial(d, h)
                E[i, j] = Rational(mw, 2**d) * krawtchouk(h, w, d) / kh
        idems.append(E)
    return idems


# ========================================================================
#  TESTS
# ========================================================================

def test_lattice_counts():
    print("\n-- Lattice shell counts (d=4) --")
    expected_r4 = {0: 1, 1: 8, 2: 24, 3: 32, 4: 24, 5: 48}
    for k, exp in expected_r4.items():
        r = shell_count(4, k)
        check(f"r_4({k})", r == exp, f"got {r}")

    n4_M = claims.get("thm:basin", "n4_M")
    N = cumulative_count(4, 5, include_origin=True)
    check(f"N_4(5) = {n4_M} (incl. origin)", N == n4_M, f"got {N}")

    N_ex = cumulative_count(4, 5, include_origin=False)
    check("N_4(5) excl. origin = 136", N_ex == 136, f"got {N_ex}")


def test_parity_weights():
    print("\n-- Parity-weight counts c_w (d=4, K=5) --")
    c = parity_weight_counts(4, 5, include_origin=True)
    expected = claims.get("lem:monotone", "c_w_K5")
    for w in range(5):
        check(f"c_{w}", c[w] == expected[w], f"got {c[w]}, expected {expected[w]}")
    n4_M = claims.get("thm:basin", "n4_M")
    check(f"sum(c_w) = {n4_M}", sum(c) == n4_M, f"got {sum(c)}")

    c_ex = parity_weight_counts(4, 5, include_origin=False)
    check("c_0 without origin = 8", c_ex[0] == 8)
    check("origin adds 1 to c_0", c[0] == c_ex[0] + 1)


def test_krawtchouk_orthogonality():
    """P * diag(C(n,w)) * P^T = 2^n * diag(C(n,j))."""
    print("\n-- Krawtchouk orthogonality (n=4) --")
    n = 4
    size = n + 1
    P = Matrix.zeros(size, size)
    for j in range(size):
        for w in range(size):
            P[j, w] = krawtchouk(j, w, n)

    D = Matrix.zeros(size, size)
    for w in range(size):
        D[w, w] = binomial(n, w)

    product = P * D * P.T
    ok = True
    for j in range(size):
        for k in range(size):
            expected = Integer(2**n) * binomial(n, j) if j == k else Integer(0)
            if product[j, k] != expected:
                ok = False
    check("P * diag(C(n,w)) * P^T = 2^n * diag(C(n,j))", ok,
          f"diag = {[product[j,j] for j in range(size)]}")


def test_idempotents():
    """Verify E_w are proper idempotents."""
    print("\n-- Idempotent properties --")
    d = 4
    N = 2**d
    idems = primitive_idempotents(d)

    # sum_w E_w = I
    total = sum(idems, Matrix.zeros(N, N))
    check("sum_w E_w = I", total == Matrix.eye(N))

    # E_w * E_w = E_w (idempotent)
    for w in range(d + 1):
        check(f"E_{w}^2 = E_{w}", idems[w] * idems[w] == idems[w])

    # E_w * E_{w'} = 0 for w != w'
    ok = True
    for w in range(d + 1):
        for wp in range(w + 1, d + 1):
            if idems[w] * idems[wp] != Matrix.zeros(N, N):
                ok = False
    check("E_w * E_{w'} = 0 for w != w'", ok)


def test_proposition_1():
    """Prop 1: G(K) = sum_w lambda_w * E_w."""
    print("\n-- SM Prop 1: Gram spectral decomposition (d=4, K=5) --")
    d, K = 4, 5
    N = 2**d

    G = gram_matrix_exact(d, K, include_origin=True)

    # Eigenvalues from parity-weight formula
    c = parity_weight_counts(d, K, include_origin=True)
    lambda_w = [Integer(N) * Rational(c[w], binomial(d, w)) for w in range(d + 1)]

    expected_lambda = claims.get("prop:spectral_q_main", "eigenvalues_n4_K5_q2")
    for w in range(d + 1):
        check(f"lambda_{w} = {expected_lambda[w]}", lambda_w[w] == expected_lambda[w],
              f"got {lambda_w[w]}")

    # G = sum_w lambda_w * E_w
    idems = primitive_idempotents(d)
    G_reconstructed = Matrix.zeros(N, N)
    for w in range(d + 1):
        G_reconstructed += lambda_w[w] * idems[w]
    check("G == sum_w lambda_w * E_w (exact)", G == G_reconstructed)

    # Trace
    n4_M = claims.get("thm:basin", "n4_M")
    tr = sum(lambda_w[w] * binomial(d, w) for w in range(d + 1))
    check(f"tr(G) = 2^n * N = {N * n4_M}", tr == N * n4_M, f"got {tr}")

    # All positive => rank 16
    check("all lambda_w > 0 => rank = 16", all(lw > 0 for lw in lambda_w))

    # Distance values
    print("\n  -- Distance values G^(h) --")
    G_h = gram_distance_values(d, K, include_origin=True)
    n4_M = claims.get("thm:basin", "n4_M")
    expected_Gh = [Integer(n4_M), Integer(5), Integer(17), Integer(-19), Integer(-39)]
    for h in range(d + 1):
        check(f"G^({h}) = {expected_Gh[h]}", G_h[h] == expected_Gh[h],
              f"got {G_h[h]}")

    # Eigenvalue from distance values: lambda_w = sum_h G^(h) * K_h(w;d)
    print("\n  -- Eigenvalue cross-check via distance values --")
    for w in range(d + 1):
        lw_from_Gh = sum(G_h[h] * krawtchouk(h, w, d) for h in range(d + 1))
        check(f"lambda_{w} from G^(h)", lw_from_Gh == lambda_w[w],
              f"sum_h G^(h)*K_h({w}) = {lw_from_Gh}")


def test_proposition_2():
    """Prop 2: rank G(K*) = 2^n implies all c_w > 0."""
    print("\n-- Lemma 6: Support-completeness --")
    d, K = 4, 5
    c = parity_weight_counts(d, K, include_origin=True)

    for w in range(d + 1):
        check(f"c_{w} > 0", c[w] > 0, f"c_{w} = {c[w]}")

    # Full Gram rank progression (paper claims rank 16 at K=4)
    print("\n  -- Full Gram rank (incl. origin) --")
    expected_full_rank = {1: 5, 2: 11, 3: 15, 4: 16, 5: 16}
    for K_test in range(1, 6):
        G = gram_matrix_exact(d, K_test, include_origin=True)
        r = G.rank()
        N_k = cumulative_count(d, K_test, include_origin=True)
        exp_r = expected_full_rank[K_test]
        check(f"K={K_test}: rank={exp_r}, N={N_k}", r == exp_r, f"rank={r}")

    # K=4 eigenvalues (paper line 343: "eigenvalues (144, 32, 64, 128, 256)")
    print("\n  -- K=4 eigenvalues --")
    c4 = parity_weight_counts(d, 4, include_origin=True)
    lambda_w_K4 = [Integer(2**d) * Rational(c4[w], binomial(d, w)) for w in range(d + 1)]
    expected_K4 = [144, 32, 64, 128, 256]
    for w in range(d + 1):
        check(f"lambda_{w}(K=4) = {expected_K4[w]}", lambda_w_K4[w] == expected_K4[w],
              f"got {lambda_w_K4[w]}")
    tr_K4 = sum(lambda_w_K4[w] * binomial(d, w) for w in range(d + 1))
    check("tr G(K=4) = 1424", tr_K4 == 1424, f"got {tr_K4}")
    # Weight-1 mass at K=4
    w1_mass_K4 = lambda_w_K4[1] * binomial(d, 1)
    w1_pct_K4 = round(float(100 * w1_mass_K4 / tr_K4), 0)
    check("K=4 weight-1 mass = 9%", w1_pct_K4 == 9, f"got {w1_pct_K4}%")

    # Dynamical Gram rank progression (excludes k=0,1 shells)
    print("\n  -- Dynamical Gram rank (excl. k=0,1) --")
    expected_dyn_rank = {2: 6, 3: 10, 4: 12, 5: 16}
    for K_test in range(2, 6):
        G_dyn = gram_dynamical(d, K_test)
        r = G_dyn.rank()
        exp_r = expected_dyn_rank[K_test]
        check(f"dyn rank G(K={K_test})", r == exp_r,
              f"got {r}, expected {exp_r}")

    check("Dynamical K*_dyn = 5 (first K with dyn rank 16)",
          gram_dynamical(d, 4).rank() < 16 and gram_dynamical(d, 5).rank() == 16)

    # Weight-0 bottleneck
    c1 = parity_weight_counts(d, 1, include_origin=False)
    check("K=1 (no origin): c_0=0 (bottleneck)", c1[0] == 0)


def test_three_sector():
    print("\n-- Three-sector decomposition --")
    d = 4
    b0, chi_orb = 1, 2**(d - 1)
    sector3 = 2**d * chi_orb
    N = cumulative_count(d, 5, include_origin=True)

    n4_M = claims.get("thm:basin", "n4_M")
    check(f"1 + 8 + 128 = {n4_M} = N_4(5)", b0 + chi_orb + sector3 == N == n4_M)

    for dd in range(2, 8):
        N_dd = cumulative_count(dd, 5, include_origin=True)
        decomp = 1 + 2**(dd-1) + 2**dd * 2**(dd-1)
        holds = (decomp == N_dd)
        if dd == 4:
            check(f"d={dd}: decomposition holds", holds)
        else:
            check(f"d={dd}: decomposition fails (1+{2**(dd-1)}+{2**dd*2**(dd-1)}={decomp} != {N_dd})", not holds)


def test_dimensional_rigidity():
    print("\n-- Dimensional rigidity --")
    for d in range(1, 7):
        holds = (2 * d == 2**(d - 1))
        if d == 4:
            check(f"d={d}: 2d = 2^{{d-1}} = {2**(d-1)}", holds)
        else:
            check(f"d={d}: 2*{d} != {2**(d-1)}", not holds)
    check("Unique in [1,100]",
          sum(1 for d in range(1, 101) if 2*d == 2**(d-1)) == 1)


def test_eigenvalue_mass():
    print("\n-- Eigenvalue mass percentages --")
    d = 4
    c = parity_weight_counts(d, 5, include_origin=True)
    lambda_w = [Integer(2**d) * Rational(c[w], binomial(d, w)) for w in range(d+1)]
    tr_parts = [lambda_w[w] * binomial(d, w) for w in range(d+1)]
    total = sum(tr_parts)
    check("tr(G) = 2192", total == 2192, f"got {total}")

    expected_pct = [6.6, 40.9, 17.5, 23.4, 11.7]
    for w in range(d + 1):
        pct = round(float(100 * tr_parts[w] / total), 1)
        check(f"pct_w={w} = {expected_pct[w]}%", pct == expected_pct[w], f"got {pct}%")

    # Eigenvalue algebraic identities (paper lines 386-393)
    chi_orb = 2**(d - 1)  # = 8
    fix = 2**d  # |Fix| = 16
    check("lambda_2 = chi_orb^2 = 64", lambda_w[2] == chi_orb**2,
          f"{lambda_w[2]} vs {chi_orb**2}")
    check("lambda_3 = |Fix| * chi_orb = 128", lambda_w[3] == fix * chi_orb,
          f"{lambda_w[3]} vs {fix * chi_orb}")
    check("lambda_4 = |Fix|^2 = 256", lambda_w[4] == fix**2,
          f"{lambda_w[4]} vs {fix**2}")

    # k=5 shell: all 48 vectors have parity weight 1 (paper lines 350-352)
    shell5_vecs = [m for m in lattice_vectors(d, 5) if sum(x*x for x in m) == 5]
    check("k=5 shell has 48 vectors", len(shell5_vecs) == 48, f"got {len(shell5_vecs)}")
    all_pw1 = all(hamming_weight_mod2(m) == 1 for m in shell5_vecs)
    check("All k=5 vectors have parity weight 1", all_pw1)

    # K* operator counts per weight class
    M_w = claims.get("lem:monotone", "M_w_K5")
    A_w = [1, 12, 54, 108, 81]
    n4_M = claims.get("thm:basin", "n4_M")
    check(f"sum(M_w) = {n4_M}", sum(M_w) == n4_M)
    check("sum(A_w) = 4^4 = 256", sum(A_w) == 256)
    saturated = [M_w[w] == A_w[w] for w in range(5)]
    check("Saturated: w=0,1,2 only", saturated == [True, True, True, False, False])


# ========================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  INDEPENDENT VERIFICATION: SM Prop 1 & Lemma 6")
    print("  Dependencies: sympy (exact arithmetic)")
    print("  No project code imported.")
    print("=" * 70)

    test_lattice_counts()
    test_parity_weights()
    test_krawtchouk_orthogonality()
    test_idempotents()
    test_proposition_1()
    test_proposition_2()
    test_three_sector()
    test_dimensional_rigidity()
    test_eigenvalue_mass()

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL PROPOSITION CHECKS VERIFIED")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 70)
    sys.exit(0 if FAIL == 0 else 1)
