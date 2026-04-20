"""Spectral Characterization of Basin Separation (Theorem 2).

Checks spectral obstruction, sharp phase transition, Delsarte
optimality, and weight saturation of K*.
"""

import numpy as np
from common import (
    pauli_tensor, pauli_weight, all_paulis, pauli_expectations,
    w_state, kstar_operator_indices, random_operator_indices,
)
from registry import claims


def _pauli_axis_state(P0_indices, n):
    """rho^+ = (I + P0)/d."""
    d = 2**n
    P0 = pauli_tensor(P0_indices)
    return (np.eye(d, dtype=complex) + P0) / d


def test_linear_bound():
    """Verify Spectral Characterization (Theorem 4) claims."""
    passed = 0
    n = 4
    d = 2**n

    # Pauli-axis states rho^pm = (I +/- P)/d have spectrum {0, 2/d}
    test_paulis = [
        (1, 0, 0, 0),
        (3, 2, 0, 0),
        (1, 2, 3, 0),
        (1, 1, 1, 1),
    ]

    for P0_idx in test_paulis:
        rho_plus = _pauli_axis_state(P0_idx, n)
        rho_minus = (np.eye(d, dtype=complex) - pauli_tensor(P0_idx)) / d

        eigs_plus = np.linalg.eigvalsh(rho_plus)
        eigs_minus = np.linalg.eigvalsh(rho_minus)

        assert np.allclose(sorted(eigs_plus), [0]*(d//2) + [2/d]*(d//2), atol=1e-12)
        assert np.all(eigs_plus >= -1e-15)
        assert np.all(eigs_minus >= -1e-15)

    passed += 1
    print(f"  [PASS] Pauli-axis states (I+/-P)/d are valid density matrices")

    # --- indistinguishability on S: test ALL unmeasured operators ---
    # Theorem 2(i): for any P0 not in S, rho+ = (I+P0)/d and rho- = (I-P0)/d
    # produce identical measurements on all P != P0 but F(rho+,rho-) = 0.

    S_kstar = kstar_operator_indices(n)
    identity = tuple([0] * n)

    # Collect all unmeasured non-identity Paulis
    unmeasured = []
    for p_idx in all_paulis(n):
        if p_idx == identity:
            continue
        if tuple(p_idx) not in S_kstar:
            unmeasured.append(tuple(p_idx))

    expected_unmeasured = claims.get("thm:spectral-char", "unmeasured_n4")
    assert len(unmeasured) == expected_unmeasured, \
        f"Expected {expected_unmeasured} unmeasured, got {len(unmeasured)}"

    # Test indistinguishability for ALL 119 unmeasured operators
    n_tested = 0
    for P0_idx in unmeasured:
        P0 = pauli_tensor(P0_idx)
        rho_plus = (np.eye(d, dtype=complex) + P0) / d
        rho_minus = (np.eye(d, dtype=complex) - P0) / d

        for Q_idx in all_paulis(n):
            if Q_idx == identity:
                continue
            Q = pauli_tensor(Q_idx)
            val_plus = np.real(np.trace(Q @ rho_plus))
            val_minus = np.real(np.trace(Q @ rho_minus))

            if tuple(Q_idx) == P0_idx:
                assert abs(val_plus - 1.0) < 1e-12
                assert abs(val_minus - (-1.0)) < 1e-12
            else:
                # Tr(Q * (I +/- P0)/d) = delta_{Q,I}/d +/- delta_{Q,P0}/d
                # For Q != I and Q != P0: both values are 0
                assert abs(val_plus) < 1e-12, \
                    f"P0={P0_idx}, Q={Q_idx}: val+ = {val_plus}"
                assert abs(val_minus) < 1e-12
        n_tested += 1

    passed += 1
    print(f"  [PASS] rho+/rho- indistinguishable on S\\{{P0}} "
          f"for ALL {n_tested} unmeasured operators")

    # F(rho+, rho-) = 0 for all unmeasured P0 (disjoint eigenvalue supports)
    n_orthogonal = 0
    for P0_idx in unmeasured:
        P0 = pauli_tensor(P0_idx)
        rho_plus = (np.eye(d, dtype=complex) + P0) / d
        rho_minus = (np.eye(d, dtype=complex) - P0) / d

        evals, evecs = np.linalg.eigh(rho_plus)
        sqrt_plus = (evecs * np.sqrt(np.maximum(evals, 0))) @ evecs.conj().T
        inner = sqrt_plus @ rho_minus @ sqrt_plus
        inner_evals = np.linalg.eigvalsh(inner)
        fidelity = np.sum(np.sqrt(np.maximum(inner_evals, 0)))**2
        assert abs(fidelity) < 1e-10, \
            f"F(rho+,rho-) = {fidelity} for P0={P0_idx}"
        n_orthogonal += 1

    passed += 1
    print(f"  [PASS] F(rho+, rho-) = 0 for all {n_orthogonal} "
          f"unmeasured operators (orthogonal)")

    # --- sharp phase transition (Theorem 2(iii)) ---
    # Compute eps_pos from actual Pauli expectations, not hardcoded values.
    # Without weight saturation: worst case F=0 (spectral obstruction above).
    # With saturation at w_sat=2: 1-F <= eps_pos = (d-1-S_k)/d^2.
    rho_W_pure = w_state(n)
    exp_W_pure = pauli_expectations(rho_W_pure, n)

    S_k = 0.0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        w = pauli_weight(p)
        x_P = exp_W_pure.get(tuple(p), 0.0)
        if w <= 2:
            S_k += x_P**2

    eps_pos_computed = (d - 1 - S_k) / d**2
    fidelity_gap = 1.0 - eps_pos_computed
    assert fidelity_gap > 0.95, \
        f"Phase transition gap = {fidelity_gap}, expected > 0.95"
    passed += 1
    print(f"  [PASS] Phase transition: 1 - eps_pos = 1 - {eps_pos_computed:.4f} "
          f"= {fidelity_gap:.4f} (computed from W-state expectations)")

    # --- Delsarte optimality ---

    import sys, os, importlib.util
    _t2_common = os.path.join(os.path.dirname(__file__), '..', 'tier2-sympy', 'common.py')
    _spec = importlib.util.spec_from_file_location("tier2_common", _t2_common)
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    qary_gram_eigenvalues = _mod.qary_gram_eigenvalues
    qary_parity_weight_counts = _mod.qary_parity_weight_counts

    # K*_spectral=4: all eigenvalues positive (support-complete)
    K_spec = 4
    eigs_kstar = qary_gram_eigenvalues(4, K_spec, 2)
    assert all(e > 0 for e in eigs_kstar), f"K*_spectral eigenvalues not all positive: {eigs_kstar}"
    passed += 1
    print(f"  [PASS] K*_spectral={K_spec}: all eigenvalues positive {[int(e) for e in eigs_kstar]}")

    # K*_spectral-1=3: zero eigenvalue (Delsarte certificate)
    eigs_below = qary_gram_eigenvalues(4, K_spec - 1, 2)
    has_zero = any(e == 0 for e in eigs_below)
    assert has_zero, f"K*_spectral-1={K_spec-1} should have zero eigenvalue: {eigs_below}"
    zero_weights = [w for w, e in enumerate(eigs_below) if e == 0]
    passed += 1
    print(f"  [PASS] K*_spectral-1={K_spec-1}: zero eigenvalue at weight(s) {zero_weights}")

    c_below = qary_parity_weight_counts(4, K_spec - 1, 2)
    for w0 in zero_weights:
        assert c_below[w0] == 0, f"c_{w0}(K*-1) = {c_below[w0]}, expected 0"
    passed += 1
    print(f"  [PASS] Delsarte certificate: c_w({K_spec-1}) = {c_below}, "
          f"zero at w={zero_weights}")

    c_kstar = qary_parity_weight_counts(4, K_spec, 2)
    assert all(c > 0 for c in c_kstar), f"c_w(K*) has zeros: {c_kstar}"
    passed += 1
    print(f"  [PASS] K*_spectral={K_spec}: all c_w > 0: {c_kstar}")

    eigs_op = qary_gram_eigenvalues(4, 5, 2)
    assert all(e > 0 for e in eigs_op)
    passed += 1
    print(f"  [PASS] K*_operational=5: all eigenvalues positive {[int(e) for e in eigs_op]}")

    # --- Delsarte universality across q ---

    # K*_spectral=4 holds universally at n=4
    for q in [3, 5]:
        eigs_at = qary_gram_eigenvalues(4, 4, q)
        eigs_below = qary_gram_eigenvalues(4, 3, q)

        assert all(e > 0 for e in eigs_at), \
            f"q={q}: K*_spectral=4 eigenvalues not all positive"
        has_zero_below = any(e == 0 for e in eigs_below)
        assert has_zero_below, \
            f"q={q}: K*_spectral-1=3 should have zero eigenvalue"

        zero_w = [w for w, e in enumerate(eigs_below) if e == 0]
        passed += 1
        print(f"  [PASS] q={q}: K*_spectral=4, Delsarte certificate at "
              f"K=3, weight(s) {zero_w}")

    # --- weight saturation of K* ---

    S_kstar = kstar_operator_indices(n)
    identity = tuple([0] * n)
    unmeasured_by_weight = {w: 0 for w in range(n + 1)}
    for P_idx in all_paulis(n):
        if P_idx == identity:
            continue
        if tuple(P_idx) not in S_kstar:
            w = pauli_weight(P_idx)
            unmeasured_by_weight[w] += 1

    assert unmeasured_by_weight[0] == 0
    assert unmeasured_by_weight[1] == 0
    assert unmeasured_by_weight[2] == 0
    assert unmeasured_by_weight[3] > 0
    assert unmeasured_by_weight[4] > 0
    passed += 1
    print(f"  [PASS] K* weight saturation: unmeasured = "
          f"{dict(unmeasured_by_weight)}")

    # P(random misses at least one weight-1 op) > 0.999
    M = len(S_kstar)
    prob_all_wt1 = 1.0
    for i in range(12):
        prob_all_wt1 *= (M - i) / (255 - i)
    assert 1 - prob_all_wt1 > 0.999
    passed += 1
    print(f"  [PASS] P(random misses wt-1) = {1-prob_all_wt1:.6f} > 0.999")

    # c_w(K) non-decreasing in K (prerequisite for Theorem 4(iv))
    for K in range(1, 6):
        c_K = qary_parity_weight_counts(4, K, 2)
        c_K1 = qary_parity_weight_counts(4, K - 1, 2)
        for w in range(5):
            assert c_K[w] >= c_K1[w], \
                f"c_{w}({K}) = {c_K[w]} < c_{w}({K-1}) = {c_K1[w]}"
    passed += 1
    print(f"  [PASS] c_w(K) monotonically non-decreasing in K")

    # --- K=4 redistribution budget (Theorem 2, Phase 2) ---
    # c_w(4) = {9,8,24,32,16}, greedy redistribution -> M_w = {1,12,28,32,16}
    from math import comb
    expected_cw4 = claims.get("lem:monotone", "c_w_K4")
    c_w_4 = qary_parity_weight_counts(4, 4, 2)
    assert list(c_w_4) == expected_cw4, \
        f"c_w(4) = {c_w_4}, expected {expected_cw4}"

    A_w = [comb(4, w) * 3**w for w in range(5)]
    # Greedy redistribution: cascade excess from low to high weight
    M_w = []
    excess = 0
    for w in range(5):
        available = c_w_4[w] + excess
        allocated = min(available, A_w[w])
        M_w.append(allocated)
        excess = available - allocated
    expected_Mw4 = claims.get("lem:monotone", "M_w_K4")
    assert M_w == expected_Mw4, \
        f"M_w(K=4) = {M_w}, expected {expected_Mw4}"
    passed += 1
    print(f"  [PASS] K=4 redistribution: c_w={list(c_w_4)}, "
          f"M_w={M_w} (52% weight-2 coverage)")

    # --- Gram condition number at operational K* ---
    # Operational K* = K*_spectral + 1 for q=2.
    # Manuscript: kappa = 4 for n=2,3,4,5; ~4.9 for n=6.
    expected_kappas = {2: 4.0, 3: 4.0, 4: 4.0, 5: 4.0, 6: 4.9}
    for n_test, expected_k in expected_kappas.items():
        # Find K*_spectral (first K with all eigenvalues > 0)
        K_spec = None
        for K_try in range(1, n_test**2 + 2):
            eigs = qary_gram_eigenvalues(n_test, K_try, 2)
            if all(e > 0 for e in eigs):
                K_spec = K_try
                break
        K_op = K_spec + 1  # operational K*
        eigs_op = qary_gram_eigenvalues(n_test, K_op, 2)
        kappa_gram = max(eigs_op) / min(eigs_op)
        assert abs(kappa_gram - expected_k) < 0.1, \
            f"n={n_test}: kappa(G(K*_op={K_op})) = {kappa_gram:.2f}, " \
            f"expected {expected_k}"
        passed += 1
        print(f"  [PASS] kappa(G(K*_op={K_op}), n={n_test}) = "
              f"{kappa_gram:.2f} (expected {expected_k})")

    # Also verify the specific n=4 claim: 256/64 = 4
    eigs_n4 = qary_gram_eigenvalues(4, 5, 2)
    expected_eigs = claims.get("prop:spectral_q_main", "eigenvalues_n4_K5_q2")
    assert max(eigs_n4) == max(expected_eigs) and min(eigs_n4) == min(expected_eigs), \
        f"n=4 K*=5 eigenvalues: max={max(eigs_n4)}, min={min(eigs_n4)}"
    passed += 1
    print(f"  [PASS] n=4 K*=5: kappa = 256/64 = 4 (manuscript Eq. kappa)")

    # --- eps_tail = 4.3e-2 for W state at k=2 (Corollary 1) ---
    rho_W = w_state(n)
    exp_W = pauli_expectations(rho_W, n)
    d_val = 2**n
    S_2 = 0.0
    for p in all_paulis(n):
        if all(x == 0 for x in p):
            continue
        w_p = pauli_weight(p)
        x_P = exp_W.get(tuple(p), 0.0)
        if w_p <= 2:
            S_2 += x_P**2
    from fractions import Fraction
    claimed_eps_tail = float(Fraction(claims.get("cor:approx_local", "eps_tail_W_k2")))
    eps_tail = (d_val - 1 - S_2) / d_val**2
    assert abs(eps_tail - claimed_eps_tail) < 1e-10, \
        f"eps_tail(W, k=2) = {eps_tail}, expected {claimed_eps_tail}"
    assert abs(eps_tail - 0.043) < 0.001, \
        f"eps_tail = {eps_tail:.4f}, manuscript claims ~4.3e-2"
    passed += 1
    print(f"  [PASS] eps_tail(W, k=2) = {eps_tail:.4f} = 11/256 "
          f"(manuscript: ~4.3e-2)")

    print(f"  [PASS] Spectral Characterization: {passed} checks")
    return passed
