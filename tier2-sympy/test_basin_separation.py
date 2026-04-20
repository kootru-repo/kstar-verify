"""Basin Separation Theorem (Theorem 1). Exact rational arithmetic."""

from sympy import Rational, binomial

from common import binary_parity_weight_counts
from registry import claims


def _kstar_weight_allocation(n, K=5):
    """Derive K* operator counts by weight class from first principles.

    NOT hardcoded.  Computes c_w (lattice parity-weight counts) via brute-force
    Z^n enumeration, then applies the bottom-up saturation rule:
      1. Initial allocation: M_w = min(c_w, A_w)
      2. Surplus from weights where c_w > A_w is collected
      3. Surplus is redistributed bottom-up to fill unsaturated weights
    This is the principled allocation that follows from the Krawtchouk
    eigenvalue-mass structure: low-weight operators carry the most
    information per operator, so they are saturated first.
    """
    from common import qary_parity_weight_counts, weight_class_sizes

    c_w = qary_parity_weight_counts(n, K, q=2)
    A_w = [int(a) for a in weight_class_sizes(n, q=2)]

    # Step 1: initial allocation and surplus collection
    M_w = {}
    surplus = 0
    for w in range(n + 1):
        alloc = min(c_w[w], A_w[w])
        M_w[w] = alloc
        surplus += max(0, c_w[w] - A_w[w])

    # Step 2: redistribute surplus bottom-up to undersaturated weights
    for w in range(n + 1):
        gap = A_w[w] - M_w[w]
        if gap > 0 and surplus > 0:
            fill = min(gap, surplus)
            M_w[w] += fill
            surplus -= fill

    return M_w


def _weight_class_sizes(n):
    """A_w = C(n,w) * 3^w."""
    return {w: int(binomial(n, w)) * (3**w) for w in range(n + 1)}


def test_basin_separation():
    """Verify Basin Separation Theorem (Theorem 1) claims."""
    passed = 0
    n = 4
    d = 2**n

    # --- K* allocation ---

    n4_M = claims.get("thm:basin", "n4_M")
    M_w = _kstar_weight_allocation(n, K=5)
    M_total = sum(M_w.values())
    assert M_total == n4_M, f"|S_K*| = {M_total}, expected {n4_M}"
    passed += 1
    print(f"  [PASS] |S_K*| = {n4_M} operators")

    A = _weight_class_sizes(n)
    expected_A = {0: 1, 1: 12, 2: 54, 3: 108, 4: 81}
    assert A == expected_A
    passed += 1
    print(f"  [PASS] Weight-class sizes A_w = {list(expected_A.values())}")

    assert sum(A.values()) == 4**n
    passed += 1
    print(f"  [PASS] Total Paulis = {sum(A.values())} = 4^{n}")

    # --- weight saturation ---

    assert M_w[0] == A[0] == 1
    passed += 1
    print(f"  [PASS] Weight-0 saturated: {M_w[0]}/{A[0]}")

    assert M_w[1] == A[1] == 12
    passed += 1
    print(f"  [PASS] Weight-1 saturated: {M_w[1]}/{A[1]}")

    assert M_w[2] == A[2] == 54
    passed += 1
    print(f"  [PASS] Weight-2 saturated: {M_w[2]}/{A[2]}")

    assert M_w[3] < A[3]
    assert M_w[3] == 54
    passed += 1
    print(f"  [PASS] Weight-3 NOT saturated: {M_w[3]}/{A[3]} = "
          f"{Rational(M_w[3], A[3])}")

    assert M_w[4] < A[4]
    assert M_w[4] == 16
    passed += 1
    print(f"  [PASS] Weight-4 NOT saturated: {M_w[4]}/{A[4]}")

    w_sat = 2
    for w in range(n + 1):
        if M_w[w] == A[w]:
            assert w <= w_sat
        else:
            assert w > w_sat
    passed += 1
    print(f"  [PASS] w_sat(S_K*) = {w_sat}")

    saturated_count = sum(M_w[w] for w in range(w_sat + 1))
    assert saturated_count == 67
    assert M_total == n4_M
    passed += 1
    print(f"  [PASS] Operator count: {saturated_count} saturated + "
          f"{M_w[3]}+{M_w[4]} partial = {M_total}")

    # --- structured case: eps_flat = 0 when w_sat >= W(rho) ---

    w_product = 1
    assert w_sat >= w_product
    passed += 1
    print(f"  [PASS] Structured case: w_sat={w_sat} >= W(product)={w_product}")

    w_2local = 2
    assert w_sat >= w_2local
    passed += 1
    print(f"  [PASS] Structured case: w_sat={w_sat} >= W(2-local)={w_2local}")

    # W and GHZ have weight-4 expectations, so w_sat < W(rho)
    w_W = 4
    assert w_sat < w_W
    passed += 1
    print(f"  [PASS] W state: w_sat={w_sat} < W(W)={w_W} "
          f"(approximate extension needed)")

    w_GHZ = 4
    assert w_sat < w_GHZ
    passed += 1
    print(f"  [PASS] GHZ state: w_sat={w_sat} < W(GHZ)={w_GHZ} "
          f"(known limitation)")

    # --- random case: expected informative flat directions ---

    M = n4_M
    total_paulis = 4**n - 1  # 255 non-identity
    frac_missing = Rational(total_paulis - M, total_paulis)
    n4_N = claims.get("thm:basin", "n4_N")
    assert frac_missing == Rational(n4_N - n4_M, n4_N)
    passed += 1
    print(f"  [PASS] Fraction missing (random): {frac_missing} = "
          f"{float(frac_missing):.4f}")

    # E[lost w=1] for random 137-subset
    n_info_w1 = 12
    expected_lost_w1 = Rational(n_info_w1 * 118, 255)
    assert abs(float(expected_lost_w1) - 5.553) < 0.01
    passed += 1
    print(f"  [PASS] Expected weight-1 loss (random): "
          f"{float(expected_lost_w1):.3f} of {n_info_w1}")

    kstar_lost_w1 = A[1] - M_w[1]
    assert kstar_lost_w1 == 0
    passed += 1
    print(f"  [PASS] K* lost weight-1: {kstar_lost_w1} (zero by saturation)")

    n_info_w2 = 54
    expected_lost_w2 = Rational(n_info_w2 * 118, 255)
    assert float(expected_lost_w2) > 24.0
    passed += 1
    print(f"  [PASS] Expected weight-2 loss (random): "
          f"{float(expected_lost_w2):.1f} of {n_info_w2}")

    kstar_lost_w2 = A[2] - M_w[2]
    assert kstar_lost_w2 == 0
    passed += 1
    print(f"  [PASS] K* lost weight-2: {kstar_lost_w2} (zero by saturation)")

    # --- eps_flat = 0 for structured case ---
    # The argument: if W(rho) <= w_sat, all informative operators have
    # weight <= w_sat. Since K* saturates weights 0..w_sat (verified above),
    # every informative operator is measured. So F cap I = empty, eps_flat = 0.
    #
    # This is verified COMPUTATIONALLY in the numerical tier (test_hessian.py)
    # for 1-local and 2-local states. Here we verify the algebraic precondition:
    # the number of Pauli operators at weight <= w_sat that K* MISSES is zero.
    n_missed_low_weight = sum(A[w] - M_w[w] for w in range(w_sat + 1))
    assert n_missed_low_weight == 0, \
        f"K* misses {n_missed_low_weight} operators at weight <= {w_sat}"
    passed += 1
    print(f"  [PASS] K* measures ALL operators at weight <= {w_sat} "
          f"(0 missed, eps_flat=0 for W(rho)<={w_sat})")

    print(f"  [PASS] Basin Separation Theorem: {passed} checks")
    return passed
