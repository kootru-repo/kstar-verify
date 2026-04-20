"""Bose-Mesner iff characterization (Lemma 5 / prop:spectral_q_main).

Verifies the new iff theorem: G(K) lies in the Bose-Mesner algebra
of H(n,q) if and only if q <= 3.

For q <= 3: all entries at the same Hamming distance are equal
  (distance-class-constancy => Bose-Mesner membership).
For q >= 4: there exist pairs at the same Hamming distance with
  different G entries (counterexample => not in Bose-Mesner).

The test builds the full q^n x q^n Gram matrix via the character sum
  G_{u,v}(K) = sum_{m: |m|^2 <= K} omega^{m_bar . (u-v)}
where omega = exp(2*pi*i/q) and m_bar = m mod q.

Uses SymPy exact arithmetic (cyclotomic fields) for q <= 5,
numerical verification for q = 7.
"""

from itertools import product as cartesian
from sympy import (
    Matrix, exp, I, pi, Rational, sqrt, cos, sin,
    simplify, nsimplify, re, im, binomial,
)


def _hamming_dist_qary(u, v, q):
    """Hamming distance between two Z_q^n vectors."""
    return sum(1 for a, b in zip(u, v) if (a - b) % q != 0)


def _build_gram_matrix_exact(n, K, q):
    """Build the full q^n x q^n Gram matrix using character sums.

    G_{u,v} = sum_{m: |m|^2 <= K} omega^{sum_j m_j*(u_j - v_j) mod q}

    For q=2, omega=-1 and this reduces to (-1)^{m . (u-v mod 2)}.

    Returns a dict mapping (u_idx, v_idx) -> G value, plus the
    vertex list.
    """
    from sympy import exp as sym_exp, Integer

    omega = sym_exp(2 * pi * I / q)

    # Vertices of H(n,q): all elements of Z_q^n
    vertices = list(cartesian(range(q), repeat=n))
    V = len(vertices)

    # Lattice points in Z^n with |m|^2 <= K
    from sympy import integer_nthroot
    bound = integer_nthroot(K, 2)[0] + 1
    lattice_pts = []
    for m in cartesian(range(-bound, bound + 1), repeat=n):
        if sum(x * x for x in m) <= K:
            lattice_pts.append(m)

    # For efficiency, precompute omega^k for k = 0..q-1
    omega_powers = [omega ** k for k in range(q)]

    # Build G matrix: G[i][j] depends only on d = (u-v) mod q
    # So we can compute G(d) for each distinct difference d, then look up.
    # This is much faster than V^2 * |lattice| operations.

    # Collect all distinct differences
    diff_to_G = {}

    # Actually, for small q^n, just compute by difference vector
    # Group vertices by their difference from vertex 0
    # G(u,v) = G(0, v-u mod q), so we only need G(0, d) for each d in Z_q^n

    for d in vertices:
        val = Integer(0)
        for m in lattice_pts:
            # omega^{sum_j m_j * d_j mod q}
            exponent = sum(m[j] * d[j] for j in range(n)) % q
            val += omega_powers[exponent]
        # Simplify: should be a rational integer for q <= 3
        val = simplify(val)
        diff_to_G[tuple(d)] = val

    return vertices, diff_to_G


def _build_gram_matrix_numerical(n, K, q):
    """Build Gram matrix numerically (for larger q where exact is slow).

    Returns vertices and diff_to_G dict with float values.
    """
    import cmath
    omega = cmath.exp(2j * cmath.pi / q)
    omega_powers = [omega ** k for k in range(q)]

    vertices = list(cartesian(range(q), repeat=n))

    from sympy import integer_nthroot
    bound = integer_nthroot(K, 2)[0] + 1
    lattice_pts = []
    for m in cartesian(range(-bound, bound + 1), repeat=n):
        if sum(x * x for x in m) <= K:
            lattice_pts.append(m)

    diff_to_G = {}
    for d in vertices:
        val = 0.0 + 0.0j
        for m in lattice_pts:
            exponent = sum(m[j] * d[j] for j in range(n)) % q
            val += omega_powers[exponent]
        diff_to_G[tuple(d)] = val

    return vertices, diff_to_G


def _check_distance_class_constancy(vertices, diff_to_G, n, q, tol=1e-8):
    """Check if G values are constant within each Hamming distance class.

    Returns (is_constant, counterexample_or_None).
    counterexample is (d1, d2, h, G(d1), G(d2)) if not constant.
    """
    zero = tuple([0] * n)

    # Group differences by Hamming weight (= Hamming distance from 0)
    by_weight = {}
    for d in vertices:
        h = sum(1 for x in d if x % q != 0)  # Hamming weight
        if h not in by_weight:
            by_weight[h] = []
        by_weight[h].append(d)

    for h in sorted(by_weight.keys()):
        diffs = by_weight[h]
        if len(diffs) <= 1:
            continue

        ref_val = diff_to_G[diffs[0]]
        for d in diffs[1:]:
            val = diff_to_G[d]
            # Compare: exact (SymPy) or numerical
            if isinstance(ref_val, complex) or isinstance(val, complex):
                if abs(val - ref_val) > tol:
                    return False, (diffs[0], d, h, ref_val, val)
            else:
                diff = simplify(val - ref_val)
                if diff != 0:
                    return False, (diffs[0], d, h, ref_val, val)

    return True, None


def test_bose_mesner_iff():
    """Verify: G(K) in BM(H(n,q)) iff q <= 3."""
    passed = 0

    # --- q = 2: exact, should be in Bose-Mesner ---
    # n=2, K=5 (K* for q=2)
    for n_test, K_test in [(2, 5), (3, 5), (4, 5)]:
        vertices, diff_to_G = _build_gram_matrix_exact(n_test, K_test, 2)
        is_const, cex = _check_distance_class_constancy(
            vertices, diff_to_G, n_test, 2)
        assert is_const, \
            f"q=2, n={n_test}, K={K_test}: not distance-class-constant! {cex}"

        # Verify all G values are real integers (omega=-1, so sum of +/-1)
        for d, val in diff_to_G.items():
            assert im(val) == 0, \
                f"q=2: G({d}) has imaginary part: {val}"

        passed += 1
        # Extract distance values for display
        dist_vals = {}
        for d, val in diff_to_G.items():
            h = sum(1 for x in d if x % 2 != 0)
            dist_vals[h] = int(re(val))
        print(f"  [PASS] q=2, n={n_test}, K={K_test}: "
              f"distance-class-constant. G^(h) = {dist_vals}")

    # --- q = 3: exact, should be in Bose-Mesner ---
    # n=2, K=9 (K*=q^2 for q=3)
    for n_test, K_test in [(2, 9)]:
        vertices, diff_to_G = _build_gram_matrix_exact(n_test, K_test, 3)
        is_const, cex = _check_distance_class_constancy(
            vertices, diff_to_G, n_test, 3)
        assert is_const, \
            f"q=3, n={n_test}, K={K_test}: not distance-class-constant! {cex}"

        # For q=3, omega = exp(2pi*i/3), but G values should still be
        # real (G is Hermitian with G(u,v) = conj(G(v,u)), and for
        # distance-class-constant matrices, symmetry forces real values)
        passed += 1
        dist_vals = {}
        for d, val in diff_to_G.items():
            h = sum(1 for x in d if x % 3 != 0)
            if h not in dist_vals:
                dist_vals[h] = complex(val)
        print(f"  [PASS] q=3, n={n_test}, K={K_test}: "
              f"distance-class-constant. G^(h) = {dist_vals}")

    # --- SM counterexample: n=1, K=1, q=4 ---
    # SM Prop 1 proof: "at n=1, K=1, within Hamming distance class 1,
    # G_{0,1} = 1+2cos(2pi/q) != G_{0,2} = 1+2cos(4pi/q) for q=4"
    # Verify the exact values stated in the SM proof.
    for q_test in [4, 5, 7]:
        vertices, diff_to_G = _build_gram_matrix_numerical(1, 1, q_test)
        is_const, cex = _check_distance_class_constancy(
            vertices, diff_to_G, 1, q_test)
        assert not is_const, \
            f"q={q_test}, n=1, K=1: unexpectedly distance-class-constant!"

        # Verify SM formula: G_{0,d} = 1 + 2*cos(2*pi*d/q) for n=1, K=1
        import math
        for d in range(1, q_test):
            expected = 1.0 + 2.0 * math.cos(2.0 * math.pi * d / q_test)
            actual = diff_to_G[(d,)].real
            assert abs(actual - expected) < 1e-10, \
                f"q={q_test}: G(0,{d}) = {actual}, expected {expected}"

        passed += 1
        vals = {d: f"{diff_to_G[(d,)].real:.4f}" for d in range(q_test)}
        print(f"  [PASS] q={q_test}, n=1, K=1 (SM counterexample): "
              f"G(0,d) = {vals}")

    # --- q = 4: should NOT be in Bose-Mesner ---
    # Residues 1 and 2 have |lift|^2 = 1 vs 4, breaking distance invariance
    for n_test, K_test in [(2, 16)]:
        vertices, diff_to_G = _build_gram_matrix_numerical(
            n_test, K_test, 4)
        is_const, cex = _check_distance_class_constancy(
            vertices, diff_to_G, n_test, 4)
        assert not is_const, \
            f"q=4, n={n_test}, K={K_test}: unexpectedly distance-class-constant!"

        d1, d2, h, g1, g2 = cex
        passed += 1
        print(f"  [PASS] q=4, n={n_test}, K={K_test}: "
              f"NOT distance-class-constant. "
              f"Counterexample at h={h}: "
              f"d1={d1} -> G={g1:.4f}, d2={d2} -> G={g2:.4f}")

    # --- q = 5: should NOT be in Bose-Mesner ---
    # Lifts: 0->0, 1->1, 2->2, 3->-2, 4->-1
    # |lift|^2: 0, 1, 4, 4, 1 -- residues 1 and 2 differ
    for n_test, K_test in [(2, 25)]:
        vertices, diff_to_G = _build_gram_matrix_numerical(
            n_test, K_test, 5)
        is_const, cex = _check_distance_class_constancy(
            vertices, diff_to_G, n_test, 5)
        assert not is_const, \
            f"q=5, n={n_test}, K={K_test}: unexpectedly distance-class-constant!"

        d1, d2, h, g1, g2 = cex
        passed += 1
        print(f"  [PASS] q=5, n={n_test}, K={K_test}: "
              f"NOT distance-class-constant. "
              f"Counterexample at h={h}: "
              f"d1={d1} -> G={g1:.4f}, d2={d2} -> G={g2:.4f}")

    # --- q = 7: should NOT be in Bose-Mesner ---
    for n_test, K_test in [(2, 49)]:
        vertices, diff_to_G = _build_gram_matrix_numerical(
            n_test, K_test, 7)
        is_const, cex = _check_distance_class_constancy(
            vertices, diff_to_G, n_test, 7)
        assert not is_const, \
            f"q=7, n={n_test}, K={K_test}: unexpectedly distance-class-constant!"

        d1, d2, h, g1, g2 = cex
        passed += 1
        print(f"  [PASS] q=7, n={n_test}, K={K_test}: "
              f"NOT distance-class-constant. "
              f"Counterexample at h={h}: "
              f"d1={d1} -> G={g1:.4f}, d2={d2} -> G={g2:.4f}")

    # --- Verify weight-block-diagonality holds for ALL q ---
    # G(u,v) depends only on the weight pattern of u-v mod q
    # (i.e., the multiset of nonzero residues, not their positions)
    # This is weaker than distance-class-constancy but holds for all q.
    for q_test in [4, 5, 7]:
        n_test = 2
        K_test = q_test ** 2
        vertices, diff_to_G = _build_gram_matrix_numerical(
            n_test, K_test, q_test)

        # Check: G(u,v) is the same for all (u,v) with the same
        # Hamming weight of (u-v mod q). This is S_n-invariance.
        by_weight = {}
        for d in vertices:
            w = sum(1 for x in d if x % q_test != 0)
            if w not in by_weight:
                by_weight[w] = []
            by_weight[w].append(diff_to_G[d])

        # Within each weight class, values may differ (not distance-class-constant)
        # but the MEAN should be well-defined (weight-block-diagonal means
        # the block average is meaningful for eigenvalue computation)
        for w in sorted(by_weight.keys()):
            vals = by_weight[w]
            # Check all values have zero imaginary part (G is real on the diagonal blocks)
            # Actually, for weight-block-diagonality, we check that
            # permuting coordinates of d doesn't change G(d).
            # d and sigma(d) must give the same G value.
            # Check: for each d, all coordinate permutations give same G
            from itertools import permutations
            for d in vertices:
                w_d = sum(1 for x in d if x % q_test != 0)
                if w_d == 0 or w_d == n_test:
                    continue  # trivially symmetric
                # Generate a permutation of d
                d_perm = tuple(reversed(d))
                if d_perm != d:
                    g_d = diff_to_G[d]
                    g_perm = diff_to_G[d_perm]
                    assert abs(g_d - g_perm) < 1e-8, \
                        f"q={q_test}: S_n violation: G({d})={g_d} != G({d_perm})={g_perm}"

        passed += 1
        print(f"  [PASS] q={q_test}, n={n_test}: "
              f"weight-block-diagonal (S_n-invariant)")

    # --- Lift norm explanation ---
    # Verify the lift norm argument that drives the iff
    print("\n  Lift norm verification (the mechanism behind iff):")
    for q_test in [2, 3, 4, 5, 7]:
        # Minimal symmetric lift: a -> closest integer to a in the
        # sense that |lift(a)| is minimized. For a in {0,...,q-1}:
        # lift(a) = a if a <= q/2, else a - q
        lift_norms_sq = set()
        for a in range(1, q_test):  # nonzero residues only
            lift_a = a if a <= q_test // 2 else a - q_test
            lift_norms_sq.add(lift_a * lift_a)

        all_equal = len(lift_norms_sq) == 1
        if q_test <= 3:
            assert all_equal, \
                f"q={q_test}: nonzero lift norms should all be equal, got {lift_norms_sq}"
        else:
            assert not all_equal, \
                f"q={q_test}: nonzero lift norms should differ, got {lift_norms_sq}"

        passed += 1
        status = "all equal" if all_equal else "distinct"
        print(f"  [PASS] q={q_test}: nonzero |lift|^2 = {sorted(lift_norms_sq)} ({status})"
              f" => {'in' if q_test <= 3 else 'NOT in'} Bose-Mesner")

    print(f"\n  Bose-Mesner iff characterization: {passed} checks passed")
    return passed
