"""Verify G(n) has full rank for n = 2, ..., 6 (gap T2-iv-2).

The Gram matrix G(n) of the K*-selected operators on the Hamming scheme
H(n, 2) has eigenvalues:

    lambda_w(n) = 2^n * c_w / C(n,w)

where c_w counts lattice points in Z^n with |m|^2 <= K* and parity weight w.

Full rank of G(n) requires all lambda_w > 0, i.e., every parity-weight class
is populated by at least one lattice point.  This is proved in the manuscript
(Prop 1) using the Krawtchouk spectral structure, but the proof is handwritten.
This test provides the direct computational check.

All arithmetic is exact (SymPy rationals).
"""

from common import (
    qary_gram_eigenvalues, binary_parity_weight_counts, weight_class_sizes,
)
from registry import claims


def test_gram_full_rank():
    """Verify all Gram eigenvalues lambda_w(n) > 0 for n = 2, ..., 6."""
    passed = 0

    # K* = n for q=2 (weight saturation at K=n)
    for n in range(2, 7):
        K = n
        eigenvalues = qary_gram_eigenvalues(n, K, q=2)
        all_positive = all(lam > 0 for lam in eigenvalues)
        assert all_positive, (
            f"G({n}) NOT full rank at K={K}: eigenvalues = {eigenvalues}"
        )
        passed += 1
        rank = sum(1 for lam in eigenvalues if lam > 0)
        print(f"  [PASS] G({n}) at K={K}: rank {rank}/{n+1}, "
              f"lambda_w = {[int(lam) for lam in eigenvalues]}")

    # Verify G(4) at K*=5 matches registry values
    ev4_k5 = qary_gram_eigenvalues(4, 5, q=2)
    c_w_k5 = binary_parity_weight_counts(4, 5)
    expected_cw = claims.get("lem:monotone", "c_w_K5")
    n4_M = claims.get("thm:basin", "n4_M")
    assert c_w_k5 == expected_cw, f"c_w(4,5) = {c_w_k5}"
    assert sum(c_w_k5) == n4_M
    passed += 1
    print(f"  [PASS] c_w(4,5) = {c_w_k5}, N_4(5) = {sum(c_w_k5)}")

    expected_eigs = claims.get("prop:spectral_q_main", "eigenvalues_n4_K5_q2")
    assert [int(e) for e in ev4_k5] == expected_eigs, (
        f"G(4) eigenvalues at K=5: {ev4_k5} != {expected_eigs}"
    )
    passed += 1
    print(f"  [PASS] G(4) at K=5: eigenvalues = {expected_eigs}")

    return passed


if __name__ == "__main__":
    total = test_gram_full_rank()
    print(f"\nTOTAL: {total} checks PASSED")
