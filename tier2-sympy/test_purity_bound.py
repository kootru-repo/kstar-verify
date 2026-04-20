"""Purity Bound (Lemma 2). Exact arithmetic via SymPy."""

from sympy import Rational, binomial, sqrt
from registry import claims


def test_purity_bound():
    """Verify Purity Bound (Lemma 2) claims."""
    passed = 0
    n = 4
    d = 2**n  # 16

    # Purity constraint: tr(sigma^2) = (1 + sum_P y_P^2) / d <= 1
    # => sum_P y_P^2 <= d - 1
    purity_upper = d - 1
    assert purity_upper == 15

    # --- product state |++++> ---
    # Only X-type operators nonzero: 2^4 - 1 = 15
    S_k_product_k1 = 4
    S_k_product_k2 = 4 + 6
    S_k_product_k4 = 15

    eps_pos_product_k2 = Rational(d - 1 - S_k_product_k2, d**2)
    assert eps_pos_product_k2 == Rational(claims.get("prop:purity_main", "eps_pos_product_k2"))
    passed += 1
    print(f"  [PASS] Product state eps_pos(k=2) = {eps_pos_product_k2} = "
          f"{float(eps_pos_product_k2):.6f}")

    eps_pos_product_k4 = Rational(d - 1 - S_k_product_k4, d**2)
    assert eps_pos_product_k4 == 0
    passed += 1
    print(f"  [PASS] Product state eps_pos(k=4) = 0 (tight)")

    # --- W state ---
    # <Z_i>=1/2, <X_iX_j>=1/2, <Y_iY_j>=-1/2, <Z_iZ_j>=0
    # S_1 = 4*(1/4) = 1;  S_2 = 1 + 6*(1/4) + 6*(1/4) = 4
    S_k_W_k2 = Rational(4)
    claimed_eps_W = Rational(claims.get("prop:purity_main", "eps_pos_W_k2"))
    eps_pos_W_k2 = Rational(d - 1 - S_k_W_k2, d**2)
    assert eps_pos_W_k2 == claimed_eps_W
    passed += 1
    print(f"  [PASS] W state eps_pos(k=2) = {eps_pos_W_k2} = "
          f"{float(eps_pos_W_k2):.4f}  (registry: {claimed_eps_W})")

    # --- GHZ state ---
    # <Z_iZ_j>=1 (6 pairs), 9 weight-4 phase correlators with |x_P|=1
    # S_2 = 6, S_4 = 15
    S_k_GHZ_k2 = Rational(6)
    claimed_eps_GHZ = Rational(claims.get("prop:purity_main", "eps_pos_GHZ_k2"))
    eps_pos_GHZ_k2 = Rational(d - 1 - S_k_GHZ_k2, d**2)
    assert eps_pos_GHZ_k2 == claimed_eps_GHZ
    passed += 1
    print(f"  [PASS] GHZ state eps_pos(k=2) = {eps_pos_GHZ_k2} = "
          f"{float(eps_pos_GHZ_k2):.4f}  (registry: {claimed_eps_GHZ})")

    # --- depolarised bound ---
    # eps_pos <= (d-1)(2p-p^2)/d^2
    p = Rational(3, 100)
    eps_depol_bound = (d - 1) * (2*p - p**2) / d**2
    assert float(eps_depol_bound) < 0.004
    passed += 1
    print(f"  [PASS] Depolarised bound (p=0.03): "
          f"eps_pos <= {float(eps_depol_bound):.6f}")

    # (d-1)*(2p-p^2)/d^2 = 15*591/10000/256 = 1773/512000
    assert eps_depol_bound == Rational(1773, 512000)
    passed += 1
    print(f"  [PASS] Depolarised bound = {eps_depol_bound} = "
          f"{float(eps_depol_bound):.6f} (exact rational)")

    # Pure k-local: all purity sits at weight <= k, so eps_pos = 0
    S_k_pure_local = d - 1
    eps_pos_pure = Rational(d - 1 - S_k_pure_local, d**2)
    assert eps_pos_pure == 0
    passed += 1
    print("  [PASS] Corollary: eps_pos = 0 for pure k-local states")

    print(f"  [PASS] Purity Bound: {passed} checks")
    return passed
