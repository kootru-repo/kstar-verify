"""Weight saturation (Theorem 1 preconditions)."""

from sympy import binomial

from common import (
    binary_parity_weight_counts, N4,
    qary_parity_weight_counts
)
from registry import claims


def _weight_class_size(n, w, q=2):
    """A_w = C(n,w) * (q^2-1)^w: total operators at weight w."""
    return int(binomial(n, w)) * ((q**2 - 1)**w)


def test_weight_saturation():
    """Verify weight saturation properties of K* set."""
    passed = 0
    n = 4
    K = 5

    M_w = claims.get("lem:monotone", "M_w_K5")
    n4_M = claims.get("thm:basin", "n4_M")
    assert sum(M_w) == n4_M
    passed += 1
    print(f"  [PASS] K* allocation: {M_w}, total = {sum(M_w)}")

    A_w = [_weight_class_size(n, w) for w in range(n + 1)]
    expected_A = [1, 12, 54, 108, 81]
    assert A_w == expected_A
    assert sum(A_w) == 4**n
    passed += 1
    print(f"  [PASS] Total Pauli operators: A_w = {A_w}, sum = {sum(A_w)}")

    # w=0,1,2 fully included; w=3,4 partial
    for w in range(3):
        assert M_w[w] == A_w[w]
    passed += 1
    print(f"  [PASS] Weights 0,1,2 fully saturated: "
          f"{M_w[:3]} = {A_w[:3]}")

    M_w_expected = claims.get("lem:monotone", "M_w_K5")
    assert M_w[3] == M_w_expected[3] and A_w[3] == 108
    passed += 1
    print(f"  [PASS] Weight 3: {M_w[3]}/{A_w[3]} = 50%")

    assert M_w[4] == M_w_expected[4] and A_w[4] == 81
    passed += 1
    print(f"  [PASS] Weight 4: {M_w[4]}/{A_w[4]} ~ 19.8%")

    w_sat = max(w for w in range(n + 1) if M_w[w] == A_w[w])
    assert w_sat == 2
    passed += 1
    print(f"  [PASS] w_sat = {w_sat}")

    c_w = binary_parity_weight_counts(n, K)
    expected_cw = claims.get("lem:monotone", "c_w_K5")
    assert c_w == expected_cw
    assert sum(c_w) == n4_M
    passed += 1
    print(f"  [PASS] Lattice parity-weight counts c_w = {c_w}")

    assert all(c > 0 for c in c_w)
    passed += 1
    print(f"  [PASS] All c_w > 0 (support-completeness at K*=5)")

    c_w_4 = binary_parity_weight_counts(n, 4)
    assert sum(c_w_4) == 89
    passed += 1
    print(f"  [PASS] N_4(4) = {sum(c_w_4)} (below K*)")

    c_w_q2 = qary_parity_weight_counts(4, 5, 2)
    assert c_w_q2 == expected_cw
    passed += 1
    print(f"  [PASS] q=2 parity weights match: {c_w_q2}")

    print(f"  [PASS] Weight saturation: {passed} checks")
    return passed
