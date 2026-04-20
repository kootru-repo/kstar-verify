"""Shared exact-arithmetic primitives for K* Verification verification.

Every function uses SymPy exact integer/rational arithmetic.
No floating point anywhere.
"""

from sympy import (
    Matrix, Rational, binomial, divisors, integer_nthroot, zeros
)
from itertools import product as cartesian


# ---- Jacobi / lattice point counting ----

def sigma_tilde_exact(k):
    """Sum of divisors of k not divisible by 4.  Exact integer."""
    if k == 0:
        return 0
    return sum(d for d in divisors(k) if d % 4 != 0)


def r4_jacobi(k):
    """r_4(k) = 8 * sigma_tilde(k) for k >= 1; r_4(0) = 1."""
    if k == 0:
        return 1
    return 8 * sigma_tilde_exact(k)


def r4_direct(k):
    """Direct brute-force count of vectors n in Z^4 with |n|^2 = k."""
    bound = integer_nthroot(k, 2)[0] + 1
    count = 0
    for n in cartesian(range(-bound, bound + 1), repeat=4):
        if sum(ni**2 for ni in n) == k:
            count += 1
    return count


def N4(K):
    """Cumulative lattice-point count: N_4(K) = sum_{k=0}^{K} r_4(k)."""
    return sum(r4_jacobi(k) for k in range(K + 1))


def Nd(d, K):
    """Cumulative lattice-point count in d dimensions.  Brute force."""
    bound = integer_nthroot(K, 2)[0] + 1
    count = 0
    for n in cartesian(range(-bound, bound + 1), repeat=d):
        if sum(ni**2 for ni in n) <= K:
            count += 1
    return count


# ---- Fixed points and lattice modes ----

def fixed_points_Z2(d):
    """Fixed points of Z_2 on T^d: all v in {0,1}^d (representing {0,1/2}^d)."""
    return list(cartesian(range(2), repeat=d))


def lattice_modes(d, K_min, K_max):
    """All n in Z^d with K_min <= |n|^2 <= K_max."""
    bound = integer_nthroot(K_max, 2)[0] + 1
    modes = []
    for n in cartesian(range(-bound, bound + 1), repeat=d):
        k = sum(ni**2 for ni in n)
        if K_min <= k <= K_max:
            modes.append(n)
    return modes


def hamming_dist(v1, v2):
    """Hamming distance between two Z_2^d vectors."""
    return sum((a - b) % 2 for a, b in zip(v1, v2))


def hamming_weight(v):
    """Hamming weight of a Z_2^d vector."""
    return sum(x % 2 for x in v)


# ---- Gram matrix ----

def gram_matrix_exact(d, K, dynamical=False):
    """Gram matrix G_{ij} = sum_n (-1)^{n.(v_i - v_j mod 2)}
    over lattice modes with |n|^2 <= K.  Exact integer.
    If dynamical=True, exclude shells k=0 and k=1.
    """
    k_min = 2 if dynamical else 0
    modes = lattice_modes(d, k_min, K)
    fps = fixed_points_Z2(d)
    n_fix = len(fps)

    G = zeros(n_fix, n_fix)
    for n in modes:
        for i in range(n_fix):
            for j in range(n_fix):
                diff = tuple((fps[i][l] - fps[j][l]) % 2 for l in range(d))
                dot_mod2 = sum(n[l] * diff[l] for l in range(d)) % 2
                G[i, j] += 1 if dot_mod2 == 0 else -1
    return G


# ---- Krawtchouk polynomials (binary) ----

def krawtchouk_exact(w, h, n=4):
    """Binary Krawtchouk polynomial K_w(h; n) over the rationals.

    K_w(h; n) = sum_{j=0}^{w} (-1)^j * C(h,j) * C(n-h, w-j)
    """
    val = Rational(0)
    for j in range(w + 1):
        val += (-1)**j * binomial(h, j) * binomial(n - h, w - j)
    return val


def krawtchouk_matrix(n=4):
    """(n+1) x (n+1) matrix [K_w(h; n)] for w,h = 0..n."""
    return Matrix(n + 1, n + 1, lambda w, h: krawtchouk_exact(w, h, n))


# ---- q-ary Krawtchouk polynomials ----

def krawtchouk_qary(k, x, n, q):
    """q-ary Krawtchouk polynomial K_k(x; n, q).

    K_k(x) = sum_{j=0}^{k} (-1)^j (q-1)^{k-j} C(x,j) C(n-x, k-j)
    """
    val = Rational(0)
    for j in range(k + 1):
        if j > x or k - j > n - x:
            continue
        val += ((-1)**j * (q - 1)**(k - j)
                * binomial(x, j) * binomial(n - x, k - j))
    return val


def krawtchouk_qary_matrix(n, q):
    """(n+1) x (n+1) q-ary Krawtchouk matrix."""
    return Matrix(n + 1, n + 1,
                  lambda w, h: krawtchouk_qary(w, h, n, q))


# ---- q-ary Gram matrix ----

def qary_parity_weight(m, q):
    """q-ary parity weight: # coordinates not divisible by q."""
    return sum(1 for x in m if x % q != 0)


def _parity_weight_slice(args):
    """Worker: count lattice points for a fixed first coordinate x0."""
    x0, n, K, q, bound = args
    residual_K = K - x0 * x0
    if residual_K < 0:
        return [0] * (n + 1)
    counts = [0] * (n + 1)
    w0 = 0 if x0 % q == 0 else 1
    for m in cartesian(range(-bound, bound + 1), repeat=n - 1):
        if sum(x * x for x in m) <= residual_K:
            w = w0 + sum(1 for x in m if x % q != 0)
            counts[w] += 1
    return counts


# Threshold: parallelize when iteration count exceeds this
_PARALLEL_THRESHOLD = 50_000_000


def qary_parity_weight_counts(n, K, q):
    """Count lattice points by q-ary parity weight.

    Returns list c_w for w = 0..n, where
    c_w = |{m in Z^n : |m|^2 <= K, pw_q(m) = w}|.

    Uses multiprocessing for large n (iteration count > 50M).
    """
    bound = integer_nthroot(K, 2)[0] + 1
    width = 2 * bound + 1
    est_iters = width ** n

    if est_iters > _PARALLEL_THRESHOLD and n >= 3:
        from multiprocessing import Pool, cpu_count
        n_workers = min(cpu_count(), width)
        args = [(x0, n, K, q, bound) for x0 in range(-bound, bound + 1)]
        with Pool(n_workers) as pool:
            results = pool.map(_parity_weight_slice, args)
        counts = [0] * (n + 1)
        for partial in results:
            for w in range(n + 1):
                counts[w] += partial[w]
        return counts

    # Small case: single-threaded
    counts = [0] * (n + 1)
    for m in cartesian(range(-bound, bound + 1), repeat=n):
        if sum(x**2 for x in m) <= K:
            w = qary_parity_weight(m, q)
            counts[w] += 1
    return counts


def qary_gram_eigenvalues(n, K, q):
    """Compute Gram eigenvalues for q-ary Hamming scheme.

    lambda_w = q^n * c_w / (C(n,w) * (q-1)^w)
    Returns list of Rational eigenvalues for w = 0..n.
    """
    c = qary_parity_weight_counts(n, K, q)
    eigs = []
    for w in range(n + 1):
        denom = binomial(n, w) * (q - 1)**w
        eigs.append(Rational(q**n * c[w], denom))
    return eigs


def qary_gram_rank_full(n, K, q):
    """Check if all q-ary Gram eigenvalues are positive (=> full rank)."""
    eigs = qary_gram_eigenvalues(n, K, q)
    return all(e > 0 for e in eigs)


def find_kstar_qary(n, q, K_max=100):
    """Find K* = min{K : all q-ary Gram eigenvalues positive}."""
    for K in range(1, K_max + 1):
        if qary_gram_rank_full(n, K, q):
            return K
    return None


# ---- Weight saturation ----

def binary_parity_weight_counts(n, K):
    """Count lattice points in Z^n with |m|^2 <= K by binary parity weight.

    Returns list c_w for w = 0..n.
    """
    return qary_parity_weight_counts(n, K, 2)


def weight_class_sizes(n, q=2):
    """Number of Pauli-type operators at each weight w = 0..n.

    For q-ary: A_w = C(n,w) * (q^2-1)^w.
    For q=2 (Pauli): A_w = C(n,w) * 3^w.
    """
    return [binomial(n, w) * (q**2 - 1)**w for w in range(n + 1)]


def weight_saturation(n, K, q=2):
    """Compute w_sat(S) for the K*-structured set on H(n,q).

    w_sat = max{k : all operators of weight <= k are included}.
    This equals max{k : sum_{w=0}^{k} A_w <= N(K)}, where
    N(K) = sum c_w = total lattice points.

    Fast path for q=2: a weight-w binary pattern is achievable iff
    there exist integers m_i with m_i odd on the w active positions
    and m_i even elsewhere, with sum m_i^2 <= K.
    Minimum cost: set active coords to +/-1 (cost 1 each), rest to 0.
    So w is saturated iff w <= K. Hence w_sat = min(K, n).
    """
    if q == 2:
        # Analytical: weight w needs min |m|^2 = w, so w_sat = min(K, n)
        return min(K, n)
    # For general q, fall back to explicit enumeration
    return _compute_wsat_explicit(n, K, q)


def _compute_wsat_explicit(n, K, q=2):
    """Compute w_sat by constructing the explicit operator set.

    For q=2: each lattice point m maps to Pauli operator with
    parity pattern (m_i mod 2). Multiple m can map to same pattern.
    Weight of the Pauli = Hamming weight of the parity pattern.
    w_sat = max{k : all weight-w Pauli patterns present for w <= k}.
    """
    bound = integer_nthroot(K, 2)[0] + 1
    # Collect all distinct parity patterns
    patterns_by_weight = {}
    for w in range(n + 1):
        patterns_by_weight[w] = set()

    for m in cartesian(range(-bound, bound + 1), repeat=n):
        if sum(x**2 for x in m) <= K:
            if q == 2:
                pattern = tuple(x % 2 for x in m)
            else:
                pattern = tuple(0 if x % q == 0 else 1 for x in m)
            w = sum(pattern)
            patterns_by_weight[w].add(pattern)

    # All C(n,w) parity patterns present => all 3^w * C(n,w) Paulis included
    w_sat = -1
    for w in range(n + 1):
        n_expected = binomial(n, w)
        if len(patterns_by_weight[w]) >= n_expected:
            w_sat = w
        else:
            break

    return w_sat
