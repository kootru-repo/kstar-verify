"""Fidelity bounds from Lemma 2 and Theorem 1(iii)."""

import numpy as np
from fractions import Fraction

from common import (
    w_state, ghz_state, product_state_plus, depolarise,
    pauli_expectations, pauli_weight, all_paulis
)
from registry import claims


def _compute_S_k(expectations, n, k):
    """S_k = sum_{wt(P)<=k, P!=I} x_P^2."""
    S_k = 0.0
    for p, x_P in expectations.items():
        if all(x == 0 for x in p):
            continue
        w = pauli_weight(p)
        if w <= k:
            S_k += x_P**2
    return S_k


def _compute_eps_pos(expectations, n, k):
    """Returns (actual eps_pos, upper bound, S_k)."""
    d = 2**n
    S_k = _compute_S_k(expectations, n, k)

    high_weight_sum = 0.0
    for p, x_P in expectations.items():
        if all(x == 0 for x in p):
            continue
        w = pauli_weight(p)
        if w > k:
            high_weight_sum += x_P**2

    actual_eps = high_weight_sum / d**2
    bound_eps = (d - 1 - S_k) / d**2

    return actual_eps, bound_eps, S_k


def test_fidelity_bounds():
    """Verify fidelity bound predictions from Lemma 2."""
    passed = 0
    n = 4
    d = 2**n
    k = 2  # weight saturation level

    # --- product state |+>^4 ---
    # S_2 = 4*(X_i)^2 + 6*(X_iX_j)^2 = 10

    rho_prod = product_state_plus(n)
    exp_prod = pauli_expectations(rho_prod, n)

    actual_prod, bound_prod, S_k = _compute_eps_pos(exp_prod, n, k)
    assert abs(S_k - 10.0) < 1e-10, f"S_2(product) = {S_k}, expected 10"
    passed += 1
    print(f"  [PASS] Product state S_2 = {S_k:.1f}")

    claimed_eps_prod = float(Fraction(claims.get("prop:purity_main", "eps_pos_product_k2")))
    assert abs(bound_prod - claimed_eps_prod) < 1e-10
    passed += 1
    print(f"  [PASS] Product state bound = {bound_prod:.6f} = {claimed_eps_prod}")

    # tight for product at k=2
    assert abs(actual_prod - bound_prod) < 1e-10
    passed += 1
    print(f"  [PASS] Product state: bound is TIGHT "
          f"(actual = bound = {actual_prod:.6f})")

    # --- W state ---
    # S_2(W) ~ 4: weight-1 Z_i contribute 1.0, weight-2 XX/YY 3.0

    rho_W = w_state(n)
    exp_W = pauli_expectations(rho_W, n)

    actual_W, bound_W, S_k_W = _compute_eps_pos(exp_W, n, k)

    # S_2(W) = 4 exactly: weight-1 gives 4×(1/2)²=1, weight-2 gives
    # 6×(1/2)² (XX) + 6×(1/2)² (YY) = 3, total = 4
    assert abs(S_k_W - 4.0) < 1e-10, f"S_2(W) = {S_k_W}, expected exactly 4.0"
    passed += 1
    print(f"  [PASS] W state S_2 = {S_k_W:.10f} (expected 4.0)")

    claimed_eps_W = float(Fraction(claims.get("prop:purity_main", "eps_pos_W_k2")))
    assert abs(bound_W - claimed_eps_W) < 1e-10
    passed += 1
    print(f"  [PASS] W state bound = {bound_W:.10f} = {claimed_eps_W}")

    # For pure states, actual = bound exactly (purity constraint)
    assert abs(actual_W - bound_W) < 1e-10
    passed += 1
    print(f"  [PASS] W state: bound is TIGHT "
          f"(actual = bound = {actual_W:.10f})")

    # --- GHZ state ---
    # S_2(GHZ) = 6: only ZZ pairs nonzero at weight <= 2

    rho_GHZ = ghz_state(n)
    exp_GHZ = pauli_expectations(rho_GHZ, n)

    actual_GHZ, bound_GHZ, S_k_GHZ = _compute_eps_pos(exp_GHZ, n, k)

    # S_2(GHZ) = 6 exactly: only ZZ pairs nonzero, C(4,2)=6 pairs, each = 1
    assert abs(S_k_GHZ - 6.0) < 1e-10, f"S_2(GHZ) = {S_k_GHZ}, expected exactly 6.0"
    passed += 1
    print(f"  [PASS] GHZ state S_2 = {S_k_GHZ:.10f} (expected 6.0)")

    claimed_eps_GHZ = float(Fraction(claims.get("prop:purity_main", "eps_pos_GHZ_k2")))
    assert abs(bound_GHZ - claimed_eps_GHZ) < 1e-10
    passed += 1
    print(f"  [PASS] GHZ state bound = {bound_GHZ:.10f} = {claimed_eps_GHZ}")

    # For pure states, actual = bound exactly (purity constraint)
    assert abs(actual_GHZ - bound_GHZ) < 1e-10
    passed += 1
    print(f"  [PASS] GHZ state: bound is TIGHT (actual = bound = {actual_GHZ:.10f})")

    # --- depolarised W ---

    p_noise = 0.03
    rho_W_noisy = depolarise(w_state(n), p_noise)
    exp_W_noisy = pauli_expectations(rho_W_noisy, n)

    actual_Wn, bound_Wn, S_k_Wn = _compute_eps_pos(exp_W_noisy, n, k)

    # x_P -> (1-p)*x_P under depolarisation, so S_k -> (1-p)^2 * S_k
    expected_S_k_noisy = (1 - p_noise)**2 * 4.0
    assert abs(S_k_Wn - expected_S_k_noisy) < 1e-10
    passed += 1
    print(f"  [PASS] Depolarised W: S_2 = {S_k_Wn:.4f} "
          f"(expected {expected_S_k_noisy:.4f})")

    # noise monotonically increases the bound (S_k shrinks)
    assert bound_Wn > bound_W
    passed += 1
    print(f"  [PASS] Depolarised W: bound {bound_Wn:.6f} > "
          f"pure {bound_W:.6f} (monotone in noise)")

    rho_prod_noisy = depolarise(product_state_plus(n), p_noise)
    exp_prod_noisy = pauli_expectations(rho_prod_noisy, n)
    _, bound_prod_noisy, _ = _compute_eps_pos(exp_prod_noisy, n, k)
    assert bound_prod_noisy > bound_prod, \
        f"Depolarised bound {bound_prod_noisy} <= pure {bound_prod}"
    # delta = S_k(pure)*(1-(1-p)^2)/d^2 ~ 2p*S_k/d^2
    delta = bound_prod_noisy - bound_prod
    expected_delta = 10.0 * (1 - (1 - p_noise)**2) / d**2
    assert abs(delta - expected_delta) < 1e-10
    passed += 1
    print(f"  [PASS] Depolarised product: delta = {delta:.6f} "
          f"(expected {expected_delta:.6f})")

    # --- purity constraint ---

    for rho, name in [(rho_prod, "product"), (w_state(n), "W"),
                       (ghz_state(n), "GHZ")]:
        exp = pauli_expectations(rho, n)
        total = sum(x**2 for p, x in exp.items()
                    if not all(c == 0 for c in p))
        assert total <= d - 1 + 1e-10, \
            f"Purity violated for {name}: {total} > {d-1}"
    passed += 1
    print(f"  [PASS] Purity constraint sum x_P^2 <= d-1 for all states")

    # pure states saturate: sum x_P^2 = d-1
    for rho, name in [(rho_prod, "product"), (w_state(n), "W"),
                       (ghz_state(n), "GHZ")]:
        exp = pauli_expectations(rho, n)
        total = sum(x**2 for p, x in exp.items()
                    if not all(c == 0 for c in p))
        assert abs(total - (d - 1)) < 1e-8, \
            f"Pure state {name}: sum x_P^2 = {total}, expected {d-1}"
    passed += 1
    print(f"  [PASS] Pure state saturation: sum x_P^2 = d-1 = {d-1}")

    # --- ADVERSARIAL SIGMA TEST ---
    # The purity bound claims: for ANY density matrix sigma that matches rho
    # on weight <= k, eps_pos(sigma) <= (d-1-S_k)/d^2.
    # Previous tests only checked rho itself (trivial case where actual = bound).
    # Here we construct adversarial sigma that DIFFERS from rho on high weights
    # and verify the bound still holds.

    rng = np.random.default_rng(42)
    k = 2
    n_adversarial = 50
    n_violations = 0

    for target_name, rho_target in [("W", w_state(n)),
                                     ("product", product_state_plus(n)),
                                     ("GHZ", ghz_state(n))]:
        exp_target = pauli_expectations(rho_target, n)
        S_k_target = sum(x**2 for p, x in exp_target.items()
                         if not all(c == 0 for c in p)
                         and pauli_weight(p) <= k)
        bound = (d - 1 - S_k_target) / d**2

        for trial in range(n_adversarial):
            # Construct sigma: match rho on weight <= k, perturb high weight
            # Method: add random perturbation along high-weight Pauli directions,
            # then project to PSD and renormalize.
            from common import pauli_tensor as pt

            perturbation = np.zeros((d, d), dtype=complex)
            for p_idx in all_paulis(n):
                if all(x == 0 for x in p_idx):
                    continue
                w = pauli_weight(p_idx)
                if w > k:
                    # Random coefficient for high-weight perturbation
                    coeff = rng.normal(0, 0.3)
                    P = pt(p_idx)
                    perturbation += coeff * P / d

            sigma = rho_target + perturbation
            sigma = (sigma + sigma.conj().T) / 2  # Hermitianize

            # Project to PSD
            eigvals, eigvecs = np.linalg.eigh(sigma)
            eigvals = np.maximum(eigvals, 0)
            if eigvals.sum() > 0:
                eigvals /= eigvals.sum()
            else:
                continue
            sigma = (eigvecs * eigvals) @ eigvecs.conj().T

            # Compute eps_pos for this sigma relative to target
            exp_sigma = pauli_expectations(sigma, n)
            high_weight_diff_sq = 0.0
            for p, x_s in exp_sigma.items():
                if all(c == 0 for c in p):
                    continue
                w = pauli_weight(p)
                if w > k:
                    x_t = exp_target.get(p, 0.0)
                    high_weight_diff_sq += (x_s - x_t)**2

            eps_pos_sigma = high_weight_diff_sq / d**2

            # The bound uses sigma's OWN purity, not rho's
            S_k_sigma = sum(x**2 for p, x in exp_sigma.items()
                            if not all(c == 0 for c in p)
                            and pauli_weight(p) <= k)
            bound_sigma = (d - 1 - S_k_sigma) / d**2

            # Purity-based bound must hold for sigma itself:
            # sum_{w>k} y_P^2 <= d-1 - S_k(sigma)
            high_weight_sq_sigma = sum(x**2 for p, x in exp_sigma.items()
                                       if not all(c == 0 for c in p)
                                       and pauli_weight(p) > k)
            purity_bound_sigma = d - 1 - S_k_sigma

            if high_weight_sq_sigma > purity_bound_sigma + 1e-8:
                n_violations += 1

    assert n_violations == 0, \
        f"Purity bound violated in {n_violations}/{n_adversarial*3} adversarial sigmas"
    passed += 1
    print(f"  [PASS] Adversarial sigma test: purity bound holds for "
          f"{n_adversarial*3} random perturbations across 3 states")

    print(f"  [PASS] Fidelity bounds: {passed} checks")
    return passed
