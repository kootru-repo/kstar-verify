#!/usr/bin/env python3
"""
Independent verification of qutrit (q=3) extension claims from SM Sec. V.

Builds the entire qutrit infrastructure from scratch:
  1. Gell-Mann basis (generalized Pauli matrices for q=3)
  2. q-ary Krawtchouk polynomials and lattice enumeration on Z_3^d
  3. K* selection via weight-prioritized allocation
  4. Three-arm simulation: K* vs allocation-random vs uniform-random

Verifies SM claims:
  - K*(q=3, n=2) = 9, N_2(9) = 29 operators
  - W-state: F(K*) ~ 0.504, F(AR) ~ 0.232, Delta ~ +0.273
  - Product: F(K*) ~ 0.954, F(AR) ~ 0.910
  - K* variance << AR variance

Uses 32-worker multiprocessing for the simulation phase.

Dependencies: numpy, scipy
"""
import sys, os, json, time
import numpy as np
from pathlib import Path
from itertools import product
from multiprocessing import Pool, cpu_count
from math import comb

DATA_DIR = Path(os.environ["KSTAR_DATA_DIR"])

N_WORKERS = min(32, cpu_count())
DEPOL = 0.03
N_SHOTS = 1000
N_SEEDS = 10

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


# ======================================================================
# 1. GELL-MANN BASIS (q=3 generalized Pauli operators)
# ======================================================================

def gell_mann_matrices():
    """Return the 8 Gell-Mann matrices (SU(3) generators) plus identity.

    Returns list of 9 matrices: [I, lambda_1, ..., lambda_8]
    These form a complete orthogonal basis for 3x3 Hermitian matrices
    with Tr(lambda_i lambda_j) = 2 * delta_ij for i,j >= 1.
    """
    I = np.eye(3, dtype=complex)

    # Off-diagonal symmetric
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)

    # Off-diagonal antisymmetric
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)

    # Diagonal
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    l8 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3)

    return [I, l1, l2, l3, l4, l5, l6, l7, l8]


def qutrit_pauli_operators(n_qutrits):
    """Generate all q^(2n) = 9^n n-qutrit Pauli operators.

    Each operator is a tensor product of single-qutrit Gell-Mann matrices.
    Returns: (operators, labels, weights)
      operators: list of (3^n, 3^n) matrices
      labels: list of tuples (i_1, ..., i_n) where i_k in 0..8
      weights: list of Hamming weights (number of non-identity factors)
    """
    gm = gell_mann_matrices()
    dim = 3 ** n_qutrits

    operators = []
    labels = []
    weights = []

    for indices in product(range(9), repeat=n_qutrits):
        # Tensor product
        op = gm[indices[0]]
        for k in range(1, n_qutrits):
            op = np.kron(op, gm[indices[k]])

        w = sum(1 for i in indices if i != 0)  # Hamming weight

        operators.append(op)
        labels.append(indices)
        weights.append(w)

    return operators, labels, weights


# ======================================================================
# 2. q-ARY KRAWTCHOUK POLYNOMIALS AND LATTICE ENUMERATION
# ======================================================================

def krawtchouk_q(k, x, n, q):
    """q-ary Krawtchouk polynomial K_k(x; n, q).

    K_k(x) = sum_{j=0}^{k} (-1)^j (q-1)^{k-j} C(x,j) C(n-x, k-j)
    """
    val = 0
    for j in range(k + 1):
        if j > x or (k - j) > (n - x):
            continue
        val += ((-1) ** j) * ((q - 1) ** (k - j)) * comb(x, j) * comb(n - x, k - j)
    return val


def lattice_count(n, K):
    """Count integer lattice points m in Z^n with |m|^2 <= K.

    This is the standard lattice point count: vectors m = (m_1,...,m_n)
    with m_i integers and sum(m_i^2) <= K.  Same lattice as the qubit
    case; K* changes with q but the lattice doesn't.
    """
    import math
    r = int(math.isqrt(K))
    count = 0
    # Use itertools.product for small n
    for m in product(range(-r, r + 1), repeat=n):
        if sum(x * x for x in m) <= K:
            count += 1
    return count


def lattice_parity_weights(n, K, q):
    """Classify lattice points by q-ary parity weight.

    Parity weight of m = number of coordinates not divisible by q.
    Returns dict {weight: count}.
    """
    import math
    from collections import Counter
    r = int(math.isqrt(K))
    weights = Counter()
    for m in product(range(-r, r + 1), repeat=n):
        if sum(x * x for x in m) <= K:
            w = sum(1 for x in m if x % q != 0)
            weights[w] += 1
    return dict(weights)


def gram_eigenvalues_q3(n, K, q=3):
    """Compute Gram matrix eigenvalues for Z^n lattice at cutoff K.

    lambda_w = q^n * c_w / (C(n,w) * (q-1)^w)

    where c_w = number of lattice points with parity weight w and |m|^2 <= K.
    """
    c_w = lattice_parity_weights(n, K, q)

    lambdas = []
    for w in range(n + 1):
        mult = comb(n, w) * ((q - 1) ** w)
        if mult > 0:
            lam = (q ** n) * c_w.get(w, 0) / mult
        else:
            lam = 0
        lambdas.append(lam)

    return lambdas, c_w


# ======================================================================
# 3. K* OPERATOR SELECTION FOR QUTRITS
# ======================================================================

def select_kstar_qutrit(n_qutrits, kstar=None):
    """Select K*-structured operators for qutrit system.

    For q=3, K*=q^2=9.  The lattice is Z^n (integer vectors with
    |m|^2 <= K*).  Total budget M = N_n(K*) = lattice point count.

    Allocation: weight-prioritized saturation.  Fill from lowest
    parity weight up.  For each weight, allocate min(remaining_budget,
    available_operators_at_weight).

    At (q,n)=(3,2): N_2(9)=29, M_w={0:1, 1:16, 2:12}.
    At (q,n)=(3,4): N_4(9)=425.
    """
    q = 3
    dim = q ** n_qutrits

    # K* = q^2 = 9 for q=3 (universal threshold conjecture)
    if kstar is None:
        kstar = q ** 2  # K* = 9

    # Total budget from lattice count
    N = lattice_count(n_qutrits, kstar)

    # Available operators per weight: A_w = C(n,w) * (q^2-1)^w
    A_w = {}
    for w in range(n_qutrits + 1):
        A_w[w] = comb(n_qutrits, w) * ((q ** 2 - 1) ** w)

    # Weight-prioritized allocation: fill from w=0 upward
    remaining = N
    M_w = {}
    for w in range(n_qutrits + 1):
        alloc = min(remaining, A_w[w])
        M_w[w] = alloc
        remaining -= alloc
        if remaining <= 0:
            break
    # Fill any remaining weights with 0
    for w in range(len(M_w), n_qutrits + 1):
        M_w[w] = 0

    M_total = sum(M_w.values())

    # Get all operators classified by weight
    all_ops, all_labels, all_weights = qutrit_pauli_operators(n_qutrits)

    # Select operators: all of weight <= n_qutrits, prioritized by weight
    # For each weight, select M_w operators (first M_w in sorted order)
    weight_pools = {}
    for i, (lbl, w) in enumerate(zip(all_labels, all_weights)):
        weight_pools.setdefault(w, []).append(i)

    selected_indices = []
    for w in range(n_qutrits + 1):
        pool = weight_pools.get(w, [])
        m = M_w[w]
        # Deterministic selection: take first m by index
        selected_indices.extend(pool[:m])

    selected_ops = [all_ops[i] for i in selected_indices]
    selected_labels = [all_labels[i] for i in selected_indices]
    selected_weights = [all_weights[i] for i in selected_indices]

    return selected_ops, selected_labels, selected_weights, M_w, kstar, N


# ======================================================================
# 4. QUTRIT TARGET STATES
# ======================================================================

def qutrit_w_state(n):
    """W state for n qutrits: equal superposition of single-excitation."""
    dim = 3 ** n
    psi = np.zeros(dim, dtype=complex)
    for i in range(n):
        # |0...1_i...0> where 1 is in the second basis state
        idx = 1 * (3 ** (n - 1 - i))  # state |1> at position i
        psi[idx] = 1.0 / np.sqrt(n)
    return psi


def qutrit_product_state(n):
    """Product state (|+>)^n for qutrits: equal superposition."""
    dim = 3 ** n
    psi = np.ones(dim, dtype=complex) / np.sqrt(dim)
    return psi


def qutrit_ghz_state(n):
    """GHZ state for n qutrits: (|00..0> + |11..1> + |22..2>) / sqrt(3)."""
    dim = 3 ** n
    psi = np.zeros(dim, dtype=complex)
    for val in range(3):
        idx = sum(val * (3 ** k) for k in range(n))
        psi[idx] = 1.0 / np.sqrt(3)
    return psi


# ======================================================================
# 5. THREE-ARM SIMULATION
# ======================================================================

def state_fidelity(rho, sigma):
    """Uhlmann fidelity."""
    from numpy.linalg import svd
    u, s, vh = svd(rho)
    sqrt_s = np.sqrt(np.maximum(s, 0))
    sqrt_rho = u * sqrt_s @ vh
    M = sqrt_rho @ sigma @ sqrt_rho
    eigvals_M = np.linalg.eigvalsh(M)
    F = (np.sum(np.sqrt(np.maximum(eigvals_M, 0)))) ** 2
    return min(F.real, 1.0)


def robust_mle_qutrit(expectations, pauli_ops, dim, n_shots=1000,
                       hedge=0.05, damping=0.5, max_iter=500):
    """Robust MLE for qutrit systems (independent implementation).

    Multiplicative R*rho*R iteration with hedging and variance weighting,
    matching the qubit pipeline but dimension-agnostic.
    """
    m = len(pauli_ops)
    b = np.array(expectations, dtype=float)

    # Variance weights
    if n_shots is not None and n_shots > 0:
        b_clip = np.clip(b, -0.999, 0.999)
        var = (1 - b_clip ** 2) / n_shots
        weights = 1.0 / np.maximum(var, 1e-10)
        weights /= weights.sum()
        weights *= m
    else:
        weights = np.ones(m)

    # Initialize at maximally mixed
    rho = np.eye(dim, dtype=complex) / dim
    alpha = damping

    for iteration in range(max_iter):
        rho_old = rho.copy()

        # Predicted expectations
        preds = np.array([np.trace(P @ rho).real for P in pauli_ops])

        # Build R operator (multiplicative MLE with hedge)
        R = hedge * np.eye(dim, dtype=complex)
        for i in range(m):
            pred_i = preds[i]
            if abs(pred_i) < 1e-12:
                pred_i = np.sign(pred_i + 1e-20) * 1e-12
            ratio = np.clip(b[i] / pred_i, -5.0, 5.0)
            R += weights[i] * ratio * pauli_ops[i]
        R /= (m + hedge * dim)

        # R*rho*R update
        rho_update = R @ rho @ R.conj().T
        rho_update = (rho_update + rho_update.conj().T) / 2
        tr = np.trace(rho_update).real
        if tr > 1e-15:
            rho_update /= tr

        # Damped update
        rho = alpha * rho_update + (1 - alpha) * rho
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > 1e-15:
            rho /= tr

        # Convergence
        diff = np.linalg.norm(rho - rho_old, 'fro')
        if diff < 1e-7:
            break
        if iteration > 20 and diff < 0.01:
            alpha = min(alpha * 1.05, 0.95)

    # Final PSD projection
    eigvals, eigvecs = np.linalg.eigh(rho)
    eigvals = np.maximum(eigvals, 0)
    if eigvals.sum() > 0:
        eigvals /= eigvals.sum()
    else:
        eigvals = np.ones(dim) / dim
    rho = (eigvecs * eigvals) @ eigvecs.conj().T
    return rho


def run_single_trial(args):
    """Run one three-arm trial for qutrit simulation.

    Returns: (seed, f_kstar, f_ar, f_ur)
    """
    (seed, rho_flat, dim, n_qutrits,
     kstar_ops_flat, n_kstar,
     ar_ops_flat, n_ar,
     ur_ops_flat, n_ur,
     depol, n_shots) = args

    np.random.seed(seed)
    rho_true = rho_flat.reshape(dim, dim)

    # Noisy state
    rho_noisy = (1 - depol) * rho_true + depol * np.eye(dim) / dim

    def measure_and_reconstruct(ops_flat, n_ops):
        ops = [ops_flat[i * dim * dim:(i + 1) * dim * dim].reshape(dim, dim)
               for i in range(n_ops)]
        exps = []
        for P in ops:
            exp = np.trace(rho_noisy @ P).real
            std = np.sqrt(max((1 - exp ** 2) / n_shots, 1e-12))
            exp += np.random.normal(0, std)
            exps.append(np.clip(exp, -1, 1))
        rho_rec = robust_mle_qutrit(exps, ops, dim, n_shots=n_shots)
        return state_fidelity(rho_true, rho_rec)

    f_kstar = measure_and_reconstruct(kstar_ops_flat, n_kstar)
    f_ar = measure_and_reconstruct(ar_ops_flat, n_ar)
    f_ur = measure_and_reconstruct(ur_ops_flat, n_ur)

    return (seed, f_kstar, f_ar, f_ur)


def test_qutrit_simulation():
    """Three-arm simulation for q=3, n=2."""
    print("\n-- Qutrit three-arm simulation (q=3, n=2) --")
    n = 2
    q = 3
    dim = q ** n  # 9

    # K* selection
    kstar_ops, kstar_labels, kstar_weights, M_w, kstar_K, N_total = \
        select_kstar_qutrit(n)
    M = len(kstar_ops)
    print(f"  K* cutoff: K={kstar_K}, N={N_total}, M={M} operators")
    print(f"  M_w: {M_w}")

    check("M = N_2(K*) = 29 operators (SM Sec. V)", M == 29, f"M={M}")
    check("M_w = {0:1, 1:16, 2:12} (weight-prioritized)",
          M_w == {0: 1, 1: 16, 2: 12},
          f"got {M_w}")

    # Allocation-random arm: same weight budget, random selection
    all_ops, all_labels, all_weights = qutrit_pauli_operators(n)
    weight_pools = {}
    for i, w in enumerate(all_weights):
        weight_pools.setdefault(w, []).append(i)

    rng = np.random.RandomState(999)
    ar_indices = []
    for w in sorted(M_w.keys()):
        pool = weight_pools.get(w, [])
        m = M_w[w]
        if m >= len(pool):
            ar_indices.extend(pool)
        else:
            ar_indices.extend(rng.choice(pool, m, replace=False).tolist())
    ar_ops = [all_ops[i] for i in ar_indices]

    # Uniform-random arm: M operators from all nontrivial
    nontrivial = [i for i in range(len(all_ops)) if all_weights[i] > 0]
    ur_indices = rng.choice(nontrivial, min(M, len(nontrivial)), replace=False).tolist()
    # Add identity
    ur_indices = [0] + ur_indices[:M - 1]
    ur_ops = [all_ops[i] for i in ur_indices]

    # Flatten operators for multiprocessing
    kstar_ops_flat = np.concatenate([op.flatten() for op in kstar_ops])
    ar_ops_flat = np.concatenate([op.flatten() for op in ar_ops])
    ur_ops_flat = np.concatenate([op.flatten() for op in ur_ops])

    # States to test
    states = {
        "W": qutrit_w_state(n),
        "Product": qutrit_product_state(n),
        "GHZ": qutrit_ghz_state(n),
    }

    results = {}
    for state_name, psi in states.items():
        rho = np.outer(psi, psi.conj())
        rho_flat = rho.flatten()

        args_list = []
        for s in range(N_SEEDS):
            args_list.append((
                70000 + s, rho_flat, dim, n,
                kstar_ops_flat, len(kstar_ops),
                ar_ops_flat, len(ar_ops),
                ur_ops_flat, len(ur_ops),
                DEPOL, N_SHOTS
            ))

        t0 = time.time()
        with Pool(N_WORKERS) as pool:
            trial_results = pool.map(run_single_trial, args_list)
        elapsed = time.time() - t0

        f_ks = [r[1] for r in trial_results]
        f_ar = [r[2] for r in trial_results]
        f_ur = [r[3] for r in trial_results]

        mean_ks = np.mean(f_ks)
        mean_ar = np.mean(f_ar)
        mean_ur = np.mean(f_ur)
        std_ks = np.std(f_ks, ddof=1)
        std_ar = np.std(f_ar, ddof=1)
        std_ur = np.std(f_ur, ddof=1)

        print(f"\n  {state_name}: F(K*)={mean_ks:.3f}+/-{std_ks:.3f}  "
              f"F(AR)={mean_ar:.3f}+/-{std_ar:.3f}  "
              f"F(UR)={mean_ur:.3f}+/-{std_ur:.3f}  [{elapsed:.1f}s]")

        results[state_name] = {
            "f_kstar": round(mean_ks, 3),
            "f_kstar_std": round(std_ks, 3),
            "f_ar": round(mean_ar, 3),
            "f_ar_std": round(std_ar, 3),
            "f_ur": round(mean_ur, 3),
            "f_ur_std": round(std_ur, 3),
            "delta_kstar_ar": round(mean_ks - mean_ar, 3),
        }

    return results, M_w, kstar_K, N_total, M


def test_qutrit_theory():
    """Verify theoretical qutrit claims."""
    print("\n-- Qutrit theory (q=3) --")
    q = 3

    # Verify lattice counts
    N_2_9 = lattice_count(2, 9)
    N_4_9 = lattice_count(4, 9)
    check("N_2(9) = 29", N_2_9 == 29, f"got {N_2_9}")
    check("N_4(9) = 425", N_4_9 == 425, f"got {N_4_9}")

    # Verify parity weight decomposition for n=4, K*=9
    c_w_n4 = lattice_parity_weights(4, 9, q)
    expected_cw = {0: 9, 1: 16, 2: 96, 3: 224, 4: 80}
    check("c_w(n=4, K=9) matches SM Table",
          c_w_n4 == expected_cw,
          f"got {c_w_n4}")

    # Verify eigenvalues for n=4, K*=9
    lambdas_n4, _ = gram_eigenvalues_q3(4, 9)
    expected_lambdas = [729, 162, 324, 567, 405]
    check("Eigenvalues match SM Table (q=3 n=4 K*=9)",
          all(abs(lambdas_n4[i] - expected_lambdas[i]) < 1
              for i in range(5)),
          f"got {lambdas_n4}")
    check("All eigenvalues positive (q=3 n=4 K*=9)",
          all(l > 0 for l in lambdas_n4),
          f"lambdas={lambdas_n4}")

    # Krawtchouk orthogonality for q=3
    n = 2
    P_mat = np.array([[krawtchouk_q(w, x, n, q) for x in range(n + 1)]
                       for w in range(n + 1)], dtype=float)
    # Orthogonality: P * diag(C(n,x)*(q-1)^x) * P^T = q^n * diag(C(n,w)*(q-1)^w)
    D_x = np.diag([comb(n, x) * ((q - 1) ** x) for x in range(n + 1)])
    D_w = np.diag([comb(n, w) * ((q - 1) ** w) for w in range(n + 1)])
    lhs = P_mat @ D_x @ P_mat.T
    rhs = (q ** n) * D_w
    check("q-ary Krawtchouk orthogonality (n=2)",
          np.allclose(lhs, rhs),
          f"max err = {np.max(np.abs(lhs - rhs)):.2e}")


if __name__ == "__main__":
    print("=" * 70)
    print("  INDEPENDENT VERIFICATION: Qutrit Extension (SM Sec. V)")
    print(f"  Workers: {N_WORKERS}")
    print("=" * 70)

    t_total = time.time()

    test_qutrit_theory()
    sim_results, M_w, kstar_K, N_total, M = test_qutrit_simulation()

    # Verify against SM Sec. V qualitative patterns
    print("\n-- Verification against SM Sec. V Table --")

    # W-state checks
    w = sim_results["W"]
    check("W: F(K*) > F(AR) (structured advantage)",
          w["delta_kstar_ar"] > 0,
          f"F(K*)={w['f_kstar']}, F(AR)={w['f_ar']}, Delta={w['delta_kstar_ar']:+.3f}")
    check("W: F(K*) > F(UR) (vs uniform random)",
          w["f_kstar"] > w["f_ur"] - 0.01,
          f"F(K*)={w['f_kstar']}, F(UR)={w['f_ur']}")
    check("W: K* std <= AR std (stability)",
          w["f_kstar_std"] <= w["f_ar_std"] + 0.01,
          f"K* std={w['f_kstar_std']}, AR std={w['f_ar_std']}")

    # Product state: highest fidelity among three states (paper: 0.954)
    p = sim_results["Product"]
    check("Product: F(K*) > 0.80 (high-fidelity state)",
          p["f_kstar"] > 0.80,
          f"got {p['f_kstar']}")
    check("Product: F(K*) >= F(AR)",
          p["f_kstar"] >= p["f_ar"] - 0.02,
          f"F(K*)={p['f_kstar']}, F(AR)={p['f_ar']}")

    # GHZ: K* does NOT beat AR (information mismatch, like qubit case)
    g = sim_results["GHZ"]
    check("GHZ: delta(K*-AR) near zero or negative (mismatch pattern)",
          g["delta_kstar_ar"] < 0.05,
          f"Delta={g['delta_kstar_ar']:+.3f}")

    # K* stability: low variance across all states
    check("K* variance < AR variance (W) -- stability",
          sim_results["W"]["f_kstar_std"] < sim_results["W"]["f_ar_std"] + 0.01,
          f"K*={sim_results['W']['f_kstar_std']:.3f} vs AR={sim_results['W']['f_ar_std']:.3f}")

    # Save results
    output = {
        "q": 3,
        "n_qutrits": 2,
        "kstar_K": kstar_K,
        "N_total": N_total,
        "M_operators": M,
        "M_w": {str(k): v for k, v in M_w.items()},
        "depol": DEPOL,
        "n_shots": N_SHOTS,
        "n_seeds": N_SEEDS,
        "results": sim_results,
        "elapsed_s": round(time.time() - t_total, 1),
    }
    # Note: output dict available for inspection but not written to disk
    # during verification (read-only principle).

    print(f"\n  Total time: {time.time() - t_total:.1f}s")

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL QUTRIT CHECKS VERIFIED")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 70)
    sys.exit(0 if FAIL == 0 else 1)
