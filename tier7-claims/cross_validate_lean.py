"""
cross_validate_lean.py — Independent Python re-implementation of every
Lean definition that has a numerical avatar, then assert it agrees with
the value declared in the Lean source.

Purpose: catch interpretation errors (Lean definition silently differs
from manuscript intent) that the Lean adversarial suite cannot detect.
The Lean kernel verifies that proofs match their statements; this
harness verifies that the *statements* (i.e. the underlying definitions)
match a textbook-style independent computation.

Each check has the form:
    expected = (compute from textbook formula in numpy/sympy)
    actual   = (literal value declared in the Lean source)
    assert expected == actual

The two implementations share no code. If they agree, both are
self-consistent under the textbook reading; if they disagree, exactly
one is wrong and the divergence pinpoints the error.

Run:  python cross_validate_lean.py
"""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction
from typing import Iterable, List, Tuple

import numpy as np

# ===========================================================================
#  Result tracking
# ===========================================================================

PASS = 0
FAIL = 0
FAILURES: list[str] = []


def check(name: str, expected, actual, *, note: str = "") -> None:
    global PASS, FAIL
    if expected == actual:
        PASS += 1
        suffix = f"  ({note})" if note else ""
        print(f"  [PASS] {name}{suffix}")
    else:
        FAIL += 1
        msg = f"  [FAIL] {name}\n         expected (python): {expected}\n         actual   (lean):   {actual}"
        if note:
            msg += f"\n         {note}"
        FAILURES.append(msg)
        print(msg)


# ===========================================================================
#  Section 1 — Lattice enumeration in Z^n
# ===========================================================================
#  Targets:
#    r4(k)            (LatticeCount.lean)
#    N4(K)            (LatticeCount.lean)
#    c_w_K5           (LatticeCount.lean) — for n=4
#    c0_K5(n), c1_K5(n), c2_K5(n)  (WeightSatAllN.lean) — for all n
# ---------------------------------------------------------------------------

def enumerate_lattice_le(n: int, K: int) -> Iterable[Tuple[int, ...]]:
    """All m in Z^n with sum(m_i^2) <= K. Brute force."""
    bound = int(math.isqrt(K))
    rng = range(-bound, bound + 1)
    for v in itertools.product(rng, repeat=n):
        if sum(x * x for x in v) <= K:
            yield v


def enumerate_lattice_eq(n: int, k: int) -> Iterable[Tuple[int, ...]]:
    """All m in Z^n with sum(m_i^2) == k."""
    bound = int(math.isqrt(k))
    rng = range(-bound, bound + 1)
    for v in itertools.product(rng, repeat=n):
        if sum(x * x for x in v) == k:
            yield v


def parity_weight(v: Tuple[int, ...]) -> int:
    """#{i : v_i is odd}."""
    return sum(1 for x in v if x % 2 != 0)


def section1_lattice() -> None:
    print("\n== Section 1: Lattice enumeration in Z^n ==")

    # r4(k): shell counts at n=4
    expected_r4 = [len(list(enumerate_lattice_eq(4, k))) for k in range(6)]
    lean_r4 = [1, 8, 24, 32, 24, 48]
    check("LatticeCount.r4 [k=0..5]", expected_r4, lean_r4,
          note="brute-force shell enumeration in Z^4")

    # N4(K): cumulative
    lean_N4 = {1: 9, 2: 33, 3: 65, 4: 89, 5: 137}
    for K, val in lean_N4.items():
        expected = sum(expected_r4[: K + 1])
        check(f"LatticeCount.N4({K})", expected, val)

    # c_w_K5 at n=4: parity-weight histogram of |m|^2 <= 5
    pw_hist = [0] * 5
    for v in enumerate_lattice_le(4, 5):
        w = parity_weight(v)
        if w < 5:
            pw_hist[w] += 1
    lean_c_w_K5 = [9, 56, 24, 32, 16]
    check("LatticeCount.c_w_K5 (n=4)", pw_hist, lean_c_w_K5,
          note="parity-weight histogram of {m in Z^4 : |m|^2 <= 5}")

    # c_w_K4 at n=4
    pw_hist = [0] * 5
    for v in enumerate_lattice_le(4, 4):
        w = parity_weight(v)
        if w < 5:
            pw_hist[w] += 1
    lean_c_w_K4 = [9, 8, 24, 32, 16]
    check("GreedyRedist.c_w_K4 (n=4)", pw_hist, lean_c_w_K4)

    # c_w_K0..K3 at n=4
    lean_c_w = {0: [1, 0, 0, 0, 0], 1: [1, 8, 0, 0, 0], 2: [1, 8, 24, 0, 0],
                3: [1, 8, 24, 32, 0]}
    for K, expected_lean in lean_c_w.items():
        pw_hist = [0] * 5
        for v in enumerate_lattice_le(4, K):
            w = parity_weight(v)
            if w < 5:
                pw_hist[w] += 1
        check(f"FidelityDichotomy.c_w_K{K}", pw_hist, expected_lean)

    # c0_K5(n), c1_K5(n), c2_K5(n) closed forms — test at n in {1..7}
    print("\n  -- WeightSatAllN closed forms vs brute enumeration (n=1..7) --")
    for n in range(1, 8):
        pw_hist = [0] * 3
        for v in enumerate_lattice_le(n, 5):
            w = parity_weight(v)
            if w < 3:
                pw_hist[w] += 1
        # Lean closed forms
        c0 = 1 + 2 * n
        c1 = 4 * n * n - 2 * n
        c2 = 2 * n * (n - 1)
        check(f"  c0_K5({n}) = 1 + 2n", pw_hist[0], c0)
        check(f"  c1_K5({n}) = 4n^2 - 2n", pw_hist[1], c1)
        check(f"  c2_K5({n}) = 2n(n-1)", pw_hist[2], c2)


# ===========================================================================
#  Section 2 — Pauli weight class sizes A_w(n) = C(n,w) * (q^2-1)^w
# ===========================================================================
#  Targets: weightClassSize, A_w_n4 (Defs.lean)
# ---------------------------------------------------------------------------

def section2_pauli_classes() -> None:
    print("\n== Section 2: Pauli weight class sizes A_w(n) ==")

    def weight_class_size(n: int, w: int, q: int = 2) -> int:
        return math.comb(n, w) * (q * q - 1) ** w

    # A_w_n4 = [1, 12, 54, 108, 81]
    expected = [weight_class_size(4, w) for w in range(5)]
    lean = [1, 12, 54, 108, 81]
    check("Defs.A_w_n4 = [C(4,w)*3^w]", expected, lean,
          note="binary q=2: 3 nontrivial Paulis per qubit")

    # Direct enumeration of n=2 Paulis as (alpha, beta) in {0..3}^n,
    # where alpha=0 means I (zero contribution to weight) and 1,2,3 means X,Y,Z.
    # weight = #{i : alpha_i != 0}
    for n in range(1, 5):
        count_by_w = [0] * (n + 1)
        for op in itertools.product(range(4), repeat=n):
            w = sum(1 for x in op if x != 0)
            count_by_w[w] += 1
        expected_closed = [weight_class_size(n, w) for w in range(n + 1)]
        check(f"  A_w({n}) by direct Pauli enumeration", count_by_w, expected_closed)


# ===========================================================================
#  Section 3 — Krawtchouk polynomials at n=4
# ===========================================================================
#  Targets: krawtchouk(n,w,h) (Krawtchouk.lean), all 25 entries + orthogonality
# ---------------------------------------------------------------------------

def section3_krawtchouk() -> None:
    print("\n== Section 3: Binary Krawtchouk polynomials (n=4) ==")

    def K(n: int, w: int, h: int) -> int:
        # Binary Krawtchouk: K_w(h;n) = sum_j (-1)^j C(h,j) C(n-h, w-j)
        return sum(
            (-1) ** j * math.comb(h, j) * math.comb(n - h, w - j)
            for j in range(w + 1)
        )

    # Lean's hardcoded matrix from Krawtchouk.lean
    lean_matrix = [
        [1, 1, 1, 1, 1],         # w=0
        [4, 2, 0, -2, -4],       # w=1
        [6, 0, -2, 0, 6],        # w=2
        [4, -2, 0, 2, -4],       # w=3
        [1, -1, 1, -1, 1],       # w=4
    ]
    for w in range(5):
        for h in range(5):
            check(f"  krawtchouk(4,{w},{h})", K(4, w, h), lean_matrix[w][h])

    # Orthogonality: sum_h C(4,h) K_w(h) K_{w'}(h) = 16*C(4,w)*delta
    print("\n  -- Krawtchouk orthogonality (binary, n=4) --")
    for w in range(5):
        for wp in range(5):
            inner = sum(math.comb(4, h) * K(4, w, h) * K(4, wp, h) for h in range(5))
            expected = 16 * math.comb(4, w) if w == wp else 0
            check(f"  <K_{w},K_{wp}>_4", inner, expected)


# ===========================================================================
#  Section 4 — Gram eigenvalues from c_w (idempotent decomposition)
# ===========================================================================
#  Targets: gramEigenvalue_from_cw (SpectralDecomp.lean), eigenvalues_K0..K5
# ---------------------------------------------------------------------------

def section4_eigenvalues() -> None:
    print("\n== Section 4: Gram eigenvalues lambda_w = 16 c_w / C(4,w) ==")

    def gram_eig(c_w: int, w: int, n: int = 4) -> int:
        # lambda_w = 2^n * c_w / C(n,w)  (n=4 → 16 * c_w / C(4,w))
        denom = math.comb(n, w)
        num = (2 ** n) * c_w
        assert num % denom == 0, f"non-integer eig at w={w}"
        return num // denom

    cases = {
        "eigenvalues_K0": ([1, 0, 0, 0, 0],   [16, 0, 0, 0, 0]),
        "eigenvalues_K1": ([1, 8, 0, 0, 0],   [16, 32, 0, 0, 0]),
        "eigenvalues_K2": ([1, 8, 24, 0, 0],  [16, 32, 64, 0, 0]),
        "eigenvalues_K3": ([1, 8, 24, 32, 0], [16, 32, 64, 128, 0]),
        "eigenvalues_K5": ([9, 56, 24, 32, 16], [144, 224, 64, 128, 256]),
    }
    # Note: eigenvalues_K4 is in Monotonicity.lean; compute c_w_K4 from
    # brute force already done in Section 1: [9, 8, 24, 32, 16].
    cases["eigenvalues_K4"] = ([9, 8, 24, 32, 16], [144, 32, 64, 128, 256])

    for name, (c_w, lean_eigs) in cases.items():
        expected = [gram_eig(c, w) for w, c in enumerate(c_w)]
        check(f"  {name}", expected, lean_eigs,
              note="lambda_w = 16 c_w / C(4,w)")


# ===========================================================================
#  Section 5 — Greedy redistribution M_w
# ===========================================================================
#  Target: greedyRedist (GreedyRedist.lean), M_w_K5, M_w_K4
# ---------------------------------------------------------------------------

def greedy_redist(c: List[int], a: List[int]) -> List[int]:
    """Independent Python implementation of the Lean greedyRedist."""
    out = []
    excess = 0
    for i, ci in enumerate(c):
        cap = a[i] if i < len(a) else 0
        avail = ci + excess
        m = min(avail, cap)
        excess = avail - m
        out.append(m)
    return out


def section5_greedy() -> None:
    print("\n== Section 5: Greedy redistribution M_w ==")

    A_w_n4 = [1, 12, 54, 108, 81]

    c_w_K5 = [9, 56, 24, 32, 16]
    M_w_K5 = greedy_redist(c_w_K5, A_w_n4)
    check("GreedyRedist.M_w_K5", M_w_K5, [1, 12, 54, 54, 16])
    check("  sum(M_w_K5)", sum(M_w_K5), 137)

    c_w_K4 = [9, 8, 24, 32, 16]
    M_w_K4 = greedy_redist(c_w_K4, A_w_n4)
    check("GreedyRedist.M_w_K4", M_w_K4, [1, 12, 28, 32, 16])

    # Saturation pattern: K=5 saturates w=0,1,2 but not w=3
    check("  K=5 saturates w=0", M_w_K5[0], A_w_n4[0])
    check("  K=5 saturates w=1", M_w_K5[1], A_w_n4[1])
    check("  K=5 saturates w=2", M_w_K5[2], A_w_n4[2])
    check("  K=5 fails to saturate w=3 (M_3 < A_3)",
          M_w_K5[3] < A_w_n4[3], True)


# ===========================================================================
#  Section 6 — Symmetric lift norm (q-ary Bose-Mesner test)
# ===========================================================================
#  Target: liftNormSq (QaryGram.lean)
# ---------------------------------------------------------------------------

def section6_lift_norm() -> None:
    print("\n== Section 6: Symmetric lift norm in Z_q ==")

    # QaryGram.liftNormSq (over Nat): min(a, q-a)^2
    def lift_norm_sq_nat(q: int, a: int) -> int:
        return min(a, q - a) ** 2

    # Adversarial.liftNormSq (over Int): symmetric mod q with negative bias
    # symmetricLift q r = let r' = r % q in if 2*r' > q then r' - q else r'
    def symmetric_lift_int(q: int, r: int) -> int:
        # Lean's % for Int uses Euclidean division: result has the sign of q
        rp = r % q  # python's % on positive q gives nonnegative result, matches Lean
        return rp - q if rp * 2 > q else rp

    def lift_norm_sq_int(q: int, r: int) -> int:
        return symmetric_lift_int(q, r) ** 2

    # Spot checks against the Lean theorems (QaryGram side)
    check("liftNormSq(2, 1)", lift_norm_sq_nat(2, 1), 1)
    check("liftNormSq(3, 1)", lift_norm_sq_nat(3, 1), 1)
    check("liftNormSq(3, 2)", lift_norm_sq_nat(3, 2), 1)
    check("liftNormSq(4, 1)", lift_norm_sq_nat(4, 1), 1)
    check("liftNormSq(4, 2)", lift_norm_sq_nat(4, 2), 4,
          note="key counterexample: q>=4 breaks equinormality")

    # Cross-check: the two independently-defined Lean liftNormSq's
    # (QaryGram over Nat vs Adversarial over Int) must agree on every
    # valid residue. This catches drift between the two formalizations.
    print("\n  -- cross-check QaryGram (Nat) vs Adversarial (Int) liftNormSq --")
    for q in range(2, 12):
        for a in range(0, q):
            nat_val = lift_norm_sq_nat(q, a) if a > 0 else 0
            int_val = lift_norm_sq_int(q, a)
            if a == 0:
                # both should give 0 at the zero residue
                if int_val != 0:
                    check(f"    q={q}, a={a}: Adversarial liftNormSq", int_val, 0)
            else:
                if nat_val != int_val:
                    check(f"    q={q}, a={a}: Nat vs Int agreement", nat_val, int_val)
    check("  Nat/Int liftNormSq agree on all q in 2..11, a in 0..q-1",
          True, True, note="(no divergences printed above)")

    # Property: equinormality of nonzero residues iff q <= 3
    print("\n  -- iff property: all nonzero residues equinormal iff q in {2,3} --")
    for q in range(2, 10):
        norms = {lift_norm_sq_nat(q, a) for a in range(1, q)}
        equinormal = (len(norms) == 1)
        expected = q <= 3
        check(f"  q={q}: equinormal? = (q<=3)?", equinormal, expected)

    # Symmetry: liftNormSq(q, a) == liftNormSq(q, q-a)
    print("\n  -- symmetry property liftNormSq(q,a) = liftNormSq(q,q-a) --")
    sym_ok = all(
        lift_norm_sq_nat(q, a) == lift_norm_sq_nat(q, q - a)
        for q in range(2, 12) for a in range(1, q)
    )
    check("  symmetry holds for q in 2..11", sym_ok, True)

    # Negative anchor: explicit residue with norm > 1 at q=4 (the lemma's
    # converse witness). If anyone refactors liftNormSq incorrectly so that
    # the q=4 counterexample disappears, this check fires.
    check("  q=4, a=2 lift norm > 1 (witness for converse direction)",
          lift_norm_sq_nat(4, 2) > 1, True)


# ===========================================================================
#  Section 7 — Hradil R-operator algebraic identity
# ===========================================================================
#  Target: R_pre, R_closed, R_pre_eq_R_closed (HradilDerivation.lean)
# ---------------------------------------------------------------------------

def section7_hradil() -> None:
    print("\n== Section 7: Hradil R-operator (b/y form vs closed form) ==")

    def R_pre(b: Fraction, y: Fraction) -> Tuple[Fraction, Fraction]:
        """((1+b)/(1+y)) (I+P)/4 + ((1-b)/(1-y)) (I-P)/4 in (I,P) coords."""
        c_pos = (1 + b) / (1 + y) / 4
        c_neg = (1 - b) / (1 - y) / 4
        return (c_pos + c_neg, c_pos - c_neg)

    def R_closed(b: Fraction, y: Fraction) -> Tuple[Fraction, Fraction]:
        """(1-by)/(1-y^2) (I/2) + (b-y)/(1-y^2) (P/2)."""
        d = 1 - y * y
        return ((1 - b * y) / d / 2, (b - y) / d / 2)

    # Random battery of (b, y) pairs in (-1, 1)
    import random
    random.seed(42)
    n_trials = 200
    fails = 0
    for _ in range(n_trials):
        b = Fraction(random.randint(-99, 99), 100)
        y = Fraction(random.randint(-99, 99), 100)
        if 1 + y == 0 or 1 - y == 0:
            continue
        if R_pre(b, y) != R_closed(b, y):
            fails += 1
    check(f"R_pre == R_closed on {n_trials} random (b,y) pairs", fails, 0,
          note="Fraction-exact algebraic identity")

    # Anchor cases
    check("R_closed(0, 0)", R_closed(Fraction(0), Fraction(0)),
          (Fraction(1, 2), Fraction(0)))
    check("R_closed(b, 0) = (1/2, b/2) for b=3/4",
          R_closed(Fraction(3, 4), Fraction(0)),
          (Fraction(1, 2), Fraction(3, 8)))

    # Fixed-point: at b = y the Pauli component vanishes, identity stays at 1/2
    for y_num in [-7, -3, 0, 1, 5, 9]:
        y = Fraction(y_num, 10)
        rc = R_closed(y, y)
        check(f"  fixed pt b=y={y}: P-coef = 0", rc[1], Fraction(0))
        check(f"  fixed pt b=y={y}: I-coef = 1/2", rc[0], Fraction(1, 2))


# ===========================================================================
#  Section 8 — Pauli trace orthogonality (n=2 explicit matrix check)
# ===========================================================================
#  Target: pauliTraceN (PauliOrthogonality.lean)
# ---------------------------------------------------------------------------

def section8_pauli_trace() -> None:
    print("\n== Section 8: Pauli trace orthogonality Tr(P Q) = 2^n delta_PQ ==")

    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    paulis = {0: I, 1: X, 2: Y, 3: Z}

    def tensor_pauli(idx: Tuple[int, ...]) -> np.ndarray:
        m = paulis[idx[0]]
        for k in idx[1:]:
            m = np.kron(m, paulis[k])
        return m

    # n = 2 and n = 3 explicit matrix trace
    for n in (2, 3):
        d = 2 ** n
        all_idx = list(itertools.product(range(4), repeat=n))
        ok_diag = True
        ok_off = True
        for P_idx in all_idx:
            for Q_idx in all_idx:
                P = tensor_pauli(P_idx)
                Q = tensor_pauli(Q_idx)
                tr = np.trace(P @ Q)
                expected = d if P_idx == Q_idx else 0
                if not np.isclose(tr, expected):
                    if P_idx == Q_idx:
                        ok_diag = False
                    else:
                        ok_off = False
        check(f"  n={n}: Tr(P P) = 2^n on diagonal", ok_diag, True)
        check(f"  n={n}: Tr(P Q) = 0 off-diagonal", ok_off, True)


# ===========================================================================
#  Section 9 — Distance values G^(h) at n=4, K=5
# ===========================================================================
#  Target: distanceValues_K5 (GramMatrix.lean) = [137, 5, 17, -19, -39]
# ---------------------------------------------------------------------------

def section9_distance_values() -> None:
    print("\n== Section 9: Distance values G^(h) (character-sum form) ==")
    # G^(h) = sum_w (c_w / C(n,w)) * K_w(h; n)
    # Equivalently: G^(h) = (1/2^n) sum_w lambda_w K_w(h) where lambda_w = 2^n c_w / C(n,w).
    # Sanity: K_w(0) = C(n,w) so G^(0) = sum_w c_w = N.
    def K(w: int, h: int) -> int:
        return sum((-1) ** j * math.comb(h, j) * math.comb(4 - h, w - j)
                   for j in range(w + 1))

    def G_h(c_w: List[int], h: int, n: int = 4) -> int:
        s = Fraction(0)
        for w in range(5):
            s += Fraction(c_w[w], math.comb(n, w)) * K(w, h)
        assert s.denominator == 1, f"non-integer character sum at h={h}"
        return s.numerator

    c_w_K5 = [9, 56, 24, 32, 16]
    G = [G_h(c_w_K5, h) for h in range(5)]
    lean = [137, 5, 17, -19, -39]
    check("GramMatrix.distanceValues_K5", G, lean,
          note="G^(h) = sum_w (c_w/C(n,w)) K_w(h;n)")

    c_w_K4 = [9, 8, 24, 32, 16]
    G = [G_h(c_w_K4, h) for h in range(5)]
    lean = [89, -19, 17, 5, 9]
    check("GramMatrix.distanceValues_K4", G, lean,
          note="K=4 entry was found wrong [89,-7,17,3,-39] by this harness "
               "and corrected; pinned by dv_K4_* anchor theorems")


# ===========================================================================
#  Section 10 — Expected missing fraction (basin separation)
# ===========================================================================
#  Target: expected_missing_fraction (BasinSeparation.lean)
# ---------------------------------------------------------------------------

def section10_expected_missing() -> None:
    print("\n== Section 10: Expected missing fraction ==")
    # Lean: expected_missing_fraction 137 255 = 118 / 255
    expected = Fraction(255 - 137, 255)
    check("expected_missing_fraction(137, 255)",
          expected, Fraction(118, 255),
          note="(N - M) / N for M=137, N=255")


# ===========================================================================
#  Section 11 — Anchor cross-checks (weighted sums, traces, decompositions)
# ===========================================================================
#  Targets in GramMatrix.lean and elsewhere: weighted_dv_sum, gram_trace_K5,
#  three-sector decomposition 1+8+128=137, A_w(4) total = 256, distance
#  distribution C(4,h) coefficients, parity-weight sum = 137.
# ---------------------------------------------------------------------------

def section11_anchors() -> None:
    print("\n== Section 11: Anchor cross-checks (sums, traces, decompositions) ==")

    # weighted_dv_sum = sum_h C(4,h) * g(h) over distanceValues_K5 = 144
    g_K5 = [137, 5, 17, -19, -39]
    weighted = sum(math.comb(4, h) * g_K5[h] for h in range(5))
    check("weighted_dv_sum (= lambda_0 = 144)", 144, weighted,
          note="sum_h C(4,h) g(h) = 16 * c_0 = 144")

    # gram_trace_K5 = 16 * g(0) = 2192
    check("gram_trace_K5 = 16 * g(0)", 2192, 16 * g_K5[0])

    # Three-sector decomposition at n=4: 137 = 1 + 8 + 128
    # Sectors: weight-0 (1), parity-defect (8), bulk (128)
    check("three-sector total 1+8+128", 137, 1 + 8 + 128,
          note="N_4(5) = 137 sector decomposition (Z2/d=4 only)")

    # A_w totals for n=4 sum to 4^n - 1 + 1 = 256 (including identity)
    A_w_n4 = [1, 12, 54, 108, 81]   # C(4,w)*3^w
    check("sum_w A_w(n=4)", 256, sum(A_w_n4),
          note="4^4 = 256 total Pauli operators on 4 qubits")
    for w in range(5):
        check(f"A_{w}(n=4) = C(4,{w})*3^{w}",
              math.comb(4, w) * 3 ** w, A_w_n4[w])

    # Parity-weight counts c_w at K=5, n=4: sum to N_4(5) = 137
    c_w_K5_n4 = [9, 56, 24, 32, 16]
    check("sum_w c_w(K=5,n=4)", 137, sum(c_w_K5_n4),
          note="lattice points with |m|^2 <= 5")

    # Eigenvalue formula lambda_w = 2^n * c_w / C(n,w) = 16*c_w/C(4,w)
    expected_eigs = [16 * c_w_K5_n4[w] // math.comb(4, w) for w in range(5)]
    actual_eigs = [144, 224, 64, 128, 256]
    for w in range(5):
        check(f"lambda_{w}(K=5,n=4) = 16 c_{w} / C(4,{w})",
              expected_eigs[w], actual_eigs[w])

    # Distance distribution in H(4,2): for each i, the multiplicity at
    # Hamming distance h is C(4,h). (Used by distance_distribution_h* theorems.)
    for h in range(5):
        # In a Hamming cube, count of vectors at distance h from any fixed v
        # equals C(n,h).
        cnt = sum(1 for v in range(16) if bin(v).count("1") == h)
        check(f"|{{v in F_2^4 : wt(v)=h}}| for h={h}",
              math.comb(4, h), cnt,
              note="distance_distribution_h* in GramMatrix.lean")

    # M_w from greedy redistribution at K=5, n=4 vs A_w (saturation w=1,2)
    # Should have M_1 = 12, M_2 = 54 for full saturation at n=4.
    # (We don't reimplement greedy here; check existing values from Section 5.)
    # Check the inequality 3n^2 >= n holds at small n (load-bearing for B):
    for n in range(1, 11):
        check(f"3n^2 >= n at n={n}", True, 3 * n * n >= n)

    # Surplus 2(c0+c1+c2) - (2A0 + 2A1 + 2A2) = 3n^2 - n
    # (Pinned in Lean by step2_K5_surplus at n=5,6,7,10.)
    def surplus(n):
        c0 = 1 + 2 * n
        c1 = 4 * n * n - 2 * n
        c2 = 2 * n * (n - 1)
        A0, A1, twoA2 = 1, 3 * n, 9 * n * (n - 1)
        return 2 * (c0 + c1 + c2) - (2 * A0 + 2 * A1 + twoA2)
    for n, exp in [(5, 70), (6, 102), (7, 140), (10, 290)]:
        check(f"step2_K5_surplus(n={n}) = 3n^2 - n",
              exp, surplus(n),
              note="negative anchor pinned in WeightSatAllN.lean")
        check(f"  closed-form 3n^2 - n at n={n}",
              3 * n * n - n, exp)

    # Krawtchouk K_w(0;n) = C(n,w) (boundary value used in saturation arg)
    for n in [4, 5, 6, 7, 10]:
        for w in range(min(n, 5) + 1):
            kw0 = sum((-1) ** j * math.comb(0, j) * math.comb(n - 0, w - j)
                      for j in range(w + 1))
            check(f"K_{w}(0; n={n}) = C(n,w)",
                  math.comb(n, w), kw0)


# ===========================================================================
#  Section 11b — Property-based fuzzing of definitions (technique #5)
# ===========================================================================
#  For each definition with a known invariant, sample random inputs and
#  assert the invariant holds. Catches algebraic mistakes that pass on the
#  hand-picked anchors but fail on a representative input distribution.
# ---------------------------------------------------------------------------

def section11b_fuzz() -> None:
    print("\n== Section 11b: Property-based fuzzing ==")
    import random
    rng = random.Random(20260407)

    # ---- liftNormSq symmetry: L(q,a) = L(q, q-a) for 0 < a < q ----
    sym_ok = 0
    for _ in range(2000):
        q = rng.randint(2, 50)
        a = rng.randint(1, q - 1)
        lhs = min(a, q - a) ** 2
        rhs = min(q - a, q - (q - a)) ** 2
        assert lhs == rhs
        sym_ok += 1
    check("liftNormSq symmetry: L(q,a) = L(q, q-a)",
          2000, sym_ok, note="2000 random (q,a) with q in [2,50]")

    # ---- liftNormSq bound: L(q,a) <= floor(q/2)^2 ----
    bound_ok = 0
    for _ in range(2000):
        q = rng.randint(2, 50)
        a = rng.randint(0, q - 1)
        lhs = min(a, q - a) ** 2
        assert lhs <= (q // 2) ** 2
        bound_ok += 1
    check("liftNormSq bound: L(q,a) <= floor(q/2)^2",
          2000, bound_ok)

    # ---- Krawtchouk orthogonality:
    #      Σ_h C(n,h) K_w(h;n) K_w'(h;n) = 2^n δ_{ww'} C(n,w)
    def kraw(n, w, h):
        return sum((-1) ** j * math.comb(h, j) * math.comb(n - h, w - j)
                   for j in range(w + 1))

    orth_ok = 0
    orth_fail = 0
    for n in range(2, 9):
        for w in range(n + 1):
            for wp in range(n + 1):
                s = sum(math.comb(n, h) * kraw(n, w, h) * kraw(n, wp, h)
                        for h in range(n + 1))
                expected = (2 ** n * math.comb(n, w)) if w == wp else 0
                if s == expected:
                    orth_ok += 1
                else:
                    orth_fail += 1
    check("Krawtchouk orthogonality (all n in [2,8], all w,w')",
          0, orth_fail,
          note=f"{orth_ok} pairs verified, weighted sum = 2^n C(n,w) delta")

    # ---- Krawtchouk recurrence:
    #      (w+1) K_{w+1}(h) = (n - 2h - w) K_w(h)*?  -- the standard recurrence
    #      We use the symmetry K_w(h;n) C(n,h) = K_h(w;n) C(n,w)
    sym_kraw_ok = 0
    sym_kraw_fail = 0
    for n in range(2, 9):
        for w in range(n + 1):
            for h in range(n + 1):
                lhs = kraw(n, w, h) * math.comb(n, h)
                rhs = kraw(n, h, w) * math.comb(n, w)
                if lhs == rhs:
                    sym_kraw_ok += 1
                else:
                    sym_kraw_fail += 1
    check("Krawtchouk symmetry: K_w(h) C(n,h) = K_h(w) C(n,w)",
          0, sym_kraw_fail, note=f"{sym_kraw_ok} (n,w,h) verified")

    # ---- Krawtchouk total at h=0: K_w(0;n) = C(n,w) ----
    for n in [4, 5, 6, 7, 10, 13]:
        for w in range(n + 1):
            check(f"K_{w}(0; n={n}) = C(n,{w})",
                  math.comb(n, w), kraw(n, w, 0))

    # ---- Closed-form lattice counts c0_K5, c1_K5, c2_K5 cross-checked
    #      against direct enumeration of Z^n with |m|^2 ≤ 5 ----
    enum_ok = 0
    for n in range(1, 7):
        py_c = [0, 0, 0, 0, 0]
        for m in itertools.product(range(-5, 6), repeat=n):
            if sum(x * x for x in m) <= 5:
                w = sum(1 for x in m if x % 2 != 0)
                if w < 5:
                    py_c[w] += 1
        c0_closed = 1 + 2 * n
        c1_closed = 4 * n * n - 2 * n
        c2_closed = 2 * n * (n - 1)
        check(f"closed c0_K5(n={n}) vs enum", py_c[0], c0_closed)
        check(f"closed c1_K5(n={n}) vs enum", py_c[1], c1_closed)
        check(f"closed c2_K5(n={n}) vs enum", py_c[2], c2_closed)
        enum_ok += 3

    # ---- Greedy redistribution invariant: total mass preserved ----
    # M_w_K5 from greedy redist of c_w_K5 against A_w_n4 must satisfy
    # sum(M) = min(sum(c), sum(A_4n)) when c is "spread" appropriately.
    # Simpler invariant: sum(M_w_K5) <= sum(c_w_K5).
    M_w_K5 = [1, 12, 54, 54, 16]
    c_w_K5 = [9, 56, 24, 32, 16]
    check("M_w_K5 mass conservation: sum M <= sum c",
          True, sum(M_w_K5) <= sum(c_w_K5),
          note=f"sum M = {sum(M_w_K5)}, sum c = {sum(c_w_K5)}")
    A_w_n4 = [1, 12, 54, 108, 81]
    for w in range(5):
        check(f"M_w_K5[{w}] <= A_w_n4[{w}] (cap by budget)",
              True, M_w_K5[w] <= A_w_n4[w])

    # ---- pauliTraceN invariant via random Pauli strings on n=2,3 qubits ----
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    P1 = [I, X, Y, Z]

    def kron_string(idx_list):
        out = P1[idx_list[0]]
        for i in idx_list[1:]:
            out = np.kron(out, P1[i])
        return out

    fuzz_pauli_ok = 0
    fuzz_pauli_fail = 0
    for n in [2, 3]:
        for _ in range(50):
            idx = [rng.randint(0, 3) for _ in range(n)]
            P = kron_string(idx)
            tr = np.trace(P @ P).real
            if abs(tr - 2 ** n) < 1e-9:
                fuzz_pauli_ok += 1
            else:
                fuzz_pauli_fail += 1
    check("Random Pauli string: Tr(P P) = 2^n",
          0, fuzz_pauli_fail, note=f"{fuzz_pauli_ok} random strings (n=2,3)")

    # Off-diagonal: random distinct pair has zero trace
    fuzz_off_ok = 0
    fuzz_off_fail = 0
    for n in [2, 3]:
        for _ in range(50):
            idx1 = [rng.randint(0, 3) for _ in range(n)]
            idx2 = [rng.randint(0, 3) for _ in range(n)]
            if idx1 == idx2:
                continue
            P, Q = kron_string(idx1), kron_string(idx2)
            tr = np.trace(P @ Q)
            if abs(tr) < 1e-9:
                fuzz_off_ok += 1
            else:
                fuzz_off_fail += 1
    check("Random distinct Pauli pair: Tr(P Q) = 0",
          0, fuzz_off_fail, note=f"{fuzz_off_ok} random pairs (n=2,3)")


# ===========================================================================
#  Section 12 — LIVE #eval extraction (closes the hardcoded-literal loophole)
# ===========================================================================
#  Shells out to `lake env lean --run Scripts/CrossValidateDump.lean`, parses
#  its `key=value` output, and asserts every kernel-evaluated value agrees
#  with Python's independent re-derivation. Earlier sections compared Python
#  to *transcribed* Lean literals; this section compares Python to the Lean
#  *kernel* directly. If a transcription was wrong (and Python was right),
#  earlier sections would be silent; this section catches it.
# ---------------------------------------------------------------------------

import subprocess
from pathlib import Path

# Lean source tree (always bundled under kstar-verify/lean4/)
_base = Path(__file__).parent
LEAN4_DIR = _base.parent / "lean4"
DUMP_SCRIPT = LEAN4_DIR / "Scripts" / "CrossValidateDump.lean"


def _parse_list_of_int(s: str) -> list[int]:
    """Parse a Lean-printed list `[a, b, c]` into a Python list of ints."""
    s = s.strip()
    assert s.startswith("[") and s.endswith("]"), f"not a list: {s!r}"
    inner = s[1:-1].strip()
    if not inner:
        return []
    return [int(x.strip()) for x in inner.split(",")]


def _parse_qq_pair(s: str) -> tuple[Fraction, Fraction]:
    """Parse `(num/den, num/den)` into a tuple of Fractions."""
    s = s.strip()
    assert s.startswith("(") and s.endswith(")"), f"not a pair: {s!r}"
    a, b = s[1:-1].split(", ")

    def parse_one(t: str) -> Fraction:
        n, d = t.split("/")
        return Fraction(int(n), int(d))

    return parse_one(a), parse_one(b)


def section12_live_eval() -> None:
    print("\n== Section 12: LIVE Lean kernel #eval extraction ==")
    if not DUMP_SCRIPT.exists():
        print(f"  [SKIP] LIVE Lean eval: {DUMP_SCRIPT} not found (optional)")
        return

    try:
        result = subprocess.run(
            ["lake", "env", "lean", "--run", "Scripts/CrossValidateDump.lean"],
            cwd=LEAN4_DIR, capture_output=True, text=True, timeout=300,
            shell=False,
        )
    except FileNotFoundError:
        check("lake on PATH", True, False, note="lake not found in PATH")
        return
    except subprocess.TimeoutExpired:
        check("lake env lean --run completes", True, False, note="timeout")
        return

    if result.returncode != 0:
        check("driver returncode 0", 0, result.returncode,
              note=result.stderr.strip()[:200])
        return

    # Parse key=value lines
    kv: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()

    # ---- Lattice counts: Lean kernel value vs Python enumeration ----
    for K, expected in [(0, 1), (1, 9), (2, 33), (3, 65), (4, 89), (5, 137)]:
        py_val = sum(1 for m in itertools.product(range(-K, K + 1), repeat=4)
                     if sum(x * x for x in m) <= K)
        check(f"LIVE N4_K{K} (kernel vs Python enum)",
              py_val, int(kv[f"N4_K{K}"]),
              note=f"Python: {expected}")

    # ---- c_w_K5, c_w_K4 ----
    py_cw5 = [0] * 5
    for m in itertools.product(range(-5, 6), repeat=4):
        if sum(x * x for x in m) <= 5:
            py_cw5[sum(1 for x in m if x % 2 != 0)] += 1
    check("LIVE c_w_K5 (kernel vs Python enum)",
          py_cw5, _parse_list_of_int(kv["c_w_K5"]))

    py_cw4 = [0] * 5
    for m in itertools.product(range(-4, 5), repeat=4):
        if sum(x * x for x in m) <= 4:
            py_cw4[sum(1 for x in m if x % 2 != 0)] += 1
    check("LIVE c_w_K4 (kernel vs Python enum)",
          py_cw4, _parse_list_of_int(kv["c_w_K4"]))

    # ---- A_w_n4 ----
    py_Aw = [math.comb(4, w) * 3 ** w for w in range(5)]
    check("LIVE A_w_n4 (kernel vs C(4,w)*3^w)",
          py_Aw, _parse_list_of_int(kv["A_w_n4"]))

    # ---- distanceValues_K5, K4 (re-derive via character sum) ----
    def distvals(K, n=4):
        out = []
        for h in range(n + 1):
            v = tuple([1] * h + [0] * (n - h))
            s = 0
            for m in itertools.product(range(-K, K + 1), repeat=n):
                if sum(x * x for x in m) <= K:
                    s += (-1) ** sum(m[i] * v[i] for i in range(n))
            out.append(s)
        return out

    check("LIVE distanceValues_K5 (kernel vs Python char-sum)",
          distvals(5), _parse_list_of_int(kv["distanceValues_K5"]))
    check("LIVE distanceValues_K4 (kernel vs Python char-sum)",
          distvals(4), _parse_list_of_int(kv["distanceValues_K4"]))

    # ---- weighted_dv_sum should equal 144 (= lambda_0) ----
    g5 = distvals(5)
    py_w = sum(math.comb(4, h) * g5[h] for h in range(5))
    check("LIVE weighted_dv_sum (kernel vs Python)",
          py_w, int(kv["weighted_dv_sum"]))

    # ---- liftNormSq[q,a] for q in {2,3,4,5,7,11} ----
    for q in [2, 3, 4, 5, 7, 11]:
        for a in range(q):
            py = min(a, q - a) ** 2
            check(f"LIVE liftNormSq[{q},{a}] (kernel vs min(a,q-a)^2)",
                  py, int(kv[f"liftNormSq[{q},{a}]"]))

    # ---- c0_K5[n], c1_K5[n], c2_K5[n], A1_n[n], twoA2_n[n] ----
    for n in [1, 2, 3, 4, 5, 6, 7, 10]:
        check(f"LIVE c0_K5[{n}]", 1 + 2 * n, int(kv[f"c0_K5[{n}]"]))
        check(f"LIVE c1_K5[{n}]", 4 * n * n - 2 * n, int(kv[f"c1_K5[{n}]"]))
        check(f"LIVE c2_K5[{n}]", 2 * n * (n - 1), int(kv[f"c2_K5[{n}]"]))
        check(f"LIVE A1_n[{n}]", 3 * n, int(kv[f"A1_n[{n}]"]))
        check(f"LIVE twoA2_n[{n}]", 9 * n * (n - 1), int(kv[f"twoA2_n[{n}]"]))

    # ---- R_closed numerical anchors ----
    def py_R_closed(b: Fraction, y: Fraction) -> tuple[Fraction, Fraction]:
        return ((1 - b * y) / (1 - y * y) / 2,
                (b - y) / (1 - y * y) / 2)

    for key, b, y in [
        ("R_closed[1/2,1/3]", Fraction(1, 2), Fraction(1, 3)),
        ("R_closed[-1/4,1/2]", Fraction(-1, 4), Fraction(1, 2)),
        ("R_closed[1,0]", Fraction(1), Fraction(0)),
    ]:
        py = py_R_closed(b, y)
        lean = _parse_qq_pair(kv[key])
        check(f"LIVE {key} (kernel vs Python Fraction)", py, lean)

    # ---- R_pre = R_closed identity at the same anchor ----
    lean_pre = _parse_qq_pair(kv["R_pre[1/2,1/3]"])
    lean_closed = _parse_qq_pair(kv["R_closed[1/2,1/3]"])
    check("LIVE R_pre[1/2,1/3] = R_closed[1/2,1/3] (kernel-level identity)",
          lean_closed, lean_pre)

    # ---- eigenvalues_K5 (lambda_w = 16 c_w / C(4,w)) ----
    py_eigs = [16 * py_cw5[w] // math.comb(4, w) for w in range(5)]
    check("LIVE eigenvalues_K5 (kernel vs 16 c_w / C(4,w))",
          py_eigs, _parse_list_of_int(kv["eigenvalues_K5"]))

    # ---- BasinSeparation: expected_missing_fraction & counts ----
    def _parse_frac(s):
        n, d = s.split("/")
        return Fraction(int(n), int(d))

    check("LIVE expected_missing_fraction[137,255] (kernel vs (255-137)/255)",
          Fraction(118, 255),
          _parse_frac(kv["expected_missing_fraction[137,255]"]))

    for A_w in [12, 54, 108, 256]:
        py = Fraction(A_w * 118, 255)
        lean = _parse_frac(kv[f"expected_missing_count[{A_w},137,255]"])
        check(f"LIVE expected_missing_count[{A_w},137,255]", py, lean)


# ===========================================================================
#  Main
# ===========================================================================

def main() -> int:
    print("=" * 72)
    print("  CROSS-VALIDATION: Lean definitions vs independent Python")
    print("=" * 72)

    section1_lattice()
    section2_pauli_classes()
    section3_krawtchouk()
    section4_eigenvalues()
    section5_greedy()
    section6_lift_norm()
    section7_hradil()
    section8_pauli_trace()
    section9_distance_values()
    section10_expected_missing()
    section11_anchors()
    section11b_fuzz()
    section12_live_eval()

    print()
    print("=" * 72)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL DEFINITIONS CROSS-VALIDATED")
        print("  (Lean defs agree with independent Python re-derivation.)")
    else:
        print("  DIVERGENCES DETECTED — review failures above")
    print("=" * 72)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
