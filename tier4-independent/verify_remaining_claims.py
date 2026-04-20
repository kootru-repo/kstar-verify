#!/usr/bin/env python3
"""
Cross-verification of remaining SM claims: ablation F values,
progressive weight-class addition, shot allocation, GHZ hardware hedge scan,
and MLE convergence behavior.

NOTE: This is CROSS-VERIFICATION, not independent verification.  It tests
the project's MLE (reconstruct_robust_mle) and noise model (depolarizing_channel)
under controlled conditions.  Independent math checks live in other scripts
(verify_propositions.py, verify_hardware_fidelities.py with cvxpy, etc.).

Imports: core.py (operator sets, states, noise), robust_mle.py (code under test)
Dependencies: numpy, multiprocessing
"""
import sys, json, time
import numpy as np
from pathlib import Path
from itertools import product as iproduct
from math import comb
from multiprocessing import Pool, cpu_count
from registry import claims

import os
# Project code under test — this script verifies MLE behavior, not math independence
from core import (select_kstar_paulis, select_random_paulis, all_pauli_operators,
                  w_state, ghz_state, state_fidelity, pauli_weight,
                  depolarizing_channel)
from robust_mle import reconstruct_robust_mle

DATA_DIR = Path(os.environ["KSTAR_DATA_DIR"])
N_WORKERS = min(32, cpu_count())
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


# ====================================================================
# Helpers
# ====================================================================

def random_pure_state_rng(n_qubits, rng):
    """Haar-random pure state as density matrix."""
    dim = 2 ** n_qubits
    psi = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def simulate_expectations(rho_noisy, ops, n_shots, rng):
    """Simulate Pauli expectations with shot noise."""
    exps = []
    for P in ops:
        true_exp = np.trace(rho_noisy @ P).real
        std = np.sqrt(max((1 - true_exp**2) / n_shots, 1e-12))
        noisy_exp = true_exp + rng.normal(0, std)
        exps.append(np.clip(noisy_exp, -1, 1))
    return exps


def get_kstar_ops_by_weight(n_qubits):
    """Get K* operators grouped by weight."""
    ops, labels, indices = select_kstar_paulis(n_qubits)
    all_ops_full, all_labels_full, all_idx_tuples = all_pauli_operators(n_qubits)
    # Map each selected operator to its weight
    weights = []
    for idx in indices:
        w = pauli_weight(all_idx_tuples[idx])
        weights.append(w)
    return ops, labels, indices, weights


# ====================================================================
# PART 1: Weight-class ablation (SM Table tab:ablation)
# ====================================================================

def _ablation_single_trial(args):
    """Single trial for ablation: generate state, test all configs."""
    state_idx, seed, ops_list, config_masks, n_shots, dep_noise = args
    rng = np.random.RandomState(seed)
    rho = random_pure_state_rng(4, rng)
    rho_noisy = depolarizing_channel(rho, dep_noise)

    results = {}
    for config_name, mask in config_masks.items():
        config_ops = [ops_list[i] for i in range(len(ops_list)) if mask[i]]
        exps = simulate_expectations(rho_noisy, config_ops, n_shots, rng)
        rho_rec = reconstruct_robust_mle(exps, config_ops, n_qubits=4,
                                         n_shots=n_shots, hedge=0.05, damping=0.5)
        f = state_fidelity(rho_rec, rho)
        results[config_name] = f
    return results


def test_ablation():
    print("\n--- PART 1: Weight-class ablation F values ---")
    t0 = time.time()

    ops, labels, indices, weights = get_kstar_ops_by_weight(4)
    n_ops = len(ops)

    # Build masks for each configuration
    configs = {}
    configs["full"] = [True] * n_ops
    for w_remove in range(5):
        name = f"remove_w{w_remove}"
        configs[name] = [weights[i] != w_remove for i in range(n_ops)]

    # Count operators per config
    for name, mask in configs.items():
        print(f"    {name}: M={sum(mask)}")

    n_states = 20
    n_noise_draws = 5
    dep_noise = 0.03
    n_shots = 1000

    # Build task list
    tasks = []
    for s in range(n_states):
        for d in range(n_noise_draws):
            seed = 42000 + s * 100 + d
            tasks.append((s, seed, ops, configs, n_shots, dep_noise))

    # Run in parallel
    with Pool(N_WORKERS) as pool:
        all_results = pool.map(_ablation_single_trial, tasks)

    # Aggregate
    f_values = {name: [] for name in configs}
    for res in all_results:
        for name, f in res.items():
            f_values[name].append(f)

    means = {name: np.mean(vals) for name, vals in f_values.items()}
    stds = {name: np.std(vals) for name, vals in f_values.items()}

    print(f"    Full K*: F = {means['full']:.3f} +/- {stds['full']:.3f}")
    for w in range(5):
        name = f"remove_w{w}"
        delta = means[name] - means['full']
        print(f"    Remove w{w}: F = {means[name]:.3f} +/- {stds[name]:.3f}, dF = {delta:+.3f}")

    # Checks
    check("Full K* F in [0.30, 0.55]",
          0.30 <= means['full'] <= 0.55,
          f"F={means['full']:.3f}")

    check("Remove w0 largest drop (dF < -0.15)",
          means['remove_w0'] - means['full'] < -0.15,
          f"dF={means['remove_w0'] - means['full']:.3f}")

    check("Remove w1 negligible (|dF| < 0.05)",
          abs(means['remove_w1'] - means['full']) < 0.05,
          f"dF={means['remove_w1'] - means['full']:.3f}")

    check("Remove w2 hurts (dF < 0)",
          means['remove_w2'] - means['full'] < 0,
          f"dF={means['remove_w2'] - means['full']:.3f}")

    check("Remove w4 neutral or helps (dF > -0.02)",
          means['remove_w4'] - means['full'] > -0.02,
          f"dF={means['remove_w4'] - means['full']:.3f}")

    print(f"    [ablation took {time.time()-t0:.1f}s]")


# ====================================================================
# PART 2: Progressive weight-class addition (SM Table tab:progressive)
# ====================================================================

def _progressive_single_trial(args):
    """Single trial for progressive: generate state, test all configs."""
    state_idx, seed, ops_list, config_masks, n_shots, dep_noise = args
    rng = np.random.RandomState(seed)
    rho = random_pure_state_rng(4, rng)
    rho_noisy = depolarizing_channel(rho, dep_noise)

    results = {}
    for config_name, mask in config_masks.items():
        config_ops = [ops_list[i] for i in range(len(ops_list)) if mask[i]]
        if len(config_ops) == 0:
            results[config_name] = 0.0
            continue
        exps = simulate_expectations(rho_noisy, config_ops, n_shots, rng)
        rho_rec = reconstruct_robust_mle(exps, config_ops, n_qubits=4,
                                         n_shots=n_shots, hedge=0.05, damping=0.5)
        f = state_fidelity(rho_rec, rho)
        results[config_name] = f
    return results


def test_progressive():
    print("\n--- PART 2: Progressive weight-class addition ---")
    t0 = time.time()

    ops, labels, indices, weights = get_kstar_ops_by_weight(4)
    n_ops = len(ops)

    # Build progressive masks: w<=0, w<=1, w<=2, w<=3, w<=4
    configs = {}
    for max_w in range(5):
        name = f"w_le_{max_w}"
        configs[name] = [weights[i] <= max_w for i in range(n_ops)]
        print(f"    {name}: M={sum(configs[name])}")

    n_states = 20
    n_noise_draws = 5
    dep_noise = 0.03
    n_shots = 1000

    tasks = []
    for s in range(n_states):
        for d in range(n_noise_draws):
            seed = 53000 + s * 100 + d
            tasks.append((s, seed, ops, configs, n_shots, dep_noise))

    with Pool(N_WORKERS) as pool:
        all_results = pool.map(_progressive_single_trial, tasks)

    f_values = {name: [] for name in configs}
    for res in all_results:
        for name, f in res.items():
            f_values[name].append(f)

    means = {name: np.mean(vals) for name, vals in f_values.items()}
    prev_f = 0.0
    deltas = {}
    for max_w in range(5):
        name = f"w_le_{max_w}"
        delta = means[name] - prev_f
        deltas[name] = delta
        print(f"    w<={max_w}: F = {means[name]:.3f}, dF = {delta:+.3f}")
        prev_f = means[name]

    # Key check: dominant jump is at w<=2 inclusion
    increments = {f"w_le_{w}": deltas[f"w_le_{w}"] for w in range(5)}
    max_jump_w = max(increments, key=increments.get)
    check("Dominant F jump at w<=2 inclusion",
          max_jump_w == "w_le_2",
          f"max jump at {max_jump_w} ({increments[max_jump_w]:+.3f})")

    # Check w<=4 adds little or reduces F
    check("w<=4 adds negligible F (|dF| < 0.05)",
          abs(deltas["w_le_4"]) < 0.05,
          f"dF={deltas['w_le_4']:+.3f}")

    print(f"    [progressive took {time.time()-t0:.1f}s]")


# ====================================================================
# PART 3: GHZ hardware hedge scan (SM Sec. XII)
# ====================================================================

def test_ghz_hedge_scan():
    print("\n--- PART 3: GHZ hardware hedge scan ---")
    t0 = time.time()

    # Load hardware data
    ghz_file = DATA_DIR / "hardware_results_ibm_fez_20260307_214922.json"
    with open(ghz_file) as f:
        data = json.load(f)

    kstar_exps = data["structured_expectations"]
    n4_M = claims.get("thm:basin", "n4_M")
    assert len(kstar_exps) == n4_M, f"Expected {n4_M}, got {len(kstar_exps)}"

    # Get K* operators
    ops, labels, indices = select_kstar_paulis(4)

    # Ideal GHZ state
    psi_ghz = ghz_state(4)
    rho_ghz = np.outer(psi_ghz, psi_ghz.conj())

    hedge_values = [0, 0.001, 0.01, 0.05, 0.1, 0.2]
    fidelities = []

    for h in hedge_values:
        rho_rec = reconstruct_robust_mle(kstar_exps, ops, n_qubits=4,
                                         n_shots=1000, hedge=max(h, 1e-6),
                                         damping=0.5)
        f = state_fidelity(rho_rec, rho_ghz)
        fidelities.append(f)
        print(f"    hedge={h:.3f}: F={f:.4f}")

    # SM claim: F ~ 0.50 for all hedge values (hedge doesn't rescue GHZ)
    # Small hedge values can cause MLE instability, so allow wider band
    all_in_range = all(0.45 <= f <= 0.55 for f in fidelities)
    check("All hedge F in [0.45, 0.55] (~ 0.50, hedge irrelevant)",
          all_in_range,
          f"range [{min(fidelities):.4f}, {max(fidelities):.4f}]")

    # Mean should be near 0.50
    f_mean = np.mean(fidelities)
    check("Mean F across hedge values near 0.50 (within 0.03)",
          abs(f_mean - 0.50) < 0.03,
          f"mean={f_mean:.4f}")

    # Key claim: no hedge value rescues fidelity above 0.55
    check("No hedge value rescues F above 0.55",
          max(fidelities) < 0.55,
          f"max F={max(fidelities):.4f}")

    print(f"    [hedge scan took {time.time()-t0:.1f}s]")


# ====================================================================
# PART 4: Shot allocation (SM Sec. IX)
# ====================================================================

def test_shot_allocation():
    print("\n--- PART 4: Shot allocation ---")
    t0 = time.time()

    # --- Arithmetic verification ---
    epsilon = 0.02
    n4_M = claims.get("thm:basin", "n4_M")
    total_budget = n4_M * 1000  # n4_M ops * 1000 shots

    # Weight-scaled allocation: shots_w proportional to 1/(1-2*eps)^(2w)
    # per operator of weight w
    ops, labels, indices, weights = get_kstar_ops_by_weight(4)
    n_ops = len(ops)

    # Compute raw allocation factors per operator
    raw_factors = np.array([1.0 / (1 - 2*epsilon)**(2*w) for w in weights])
    # Normalize so total shots = total_budget
    norm_factor = total_budget / raw_factors.sum()
    shots_per_op = raw_factors * norm_factor

    # Expected per-weight-class shots (per operator)
    weight_classes = sorted(set(weights))
    print("    Weight-scaled shots per operator:")
    for w in weight_classes:
        mask = [i for i, ww in enumerate(weights) if ww == w]
        avg_shots = np.mean(shots_per_op[mask])
        print(f"      w={w}: {avg_shots:.0f} shots/op (n={len(mask)} ops)")

    # SM claims approximately: 812 (w0) to 1124 (w4), ratio ~1.38x
    w0_shots = np.mean([shots_per_op[i] for i in range(n_ops) if weights[i] == 0])
    w4_shots = np.mean([shots_per_op[i] for i in range(n_ops) if weights[i] == 4])
    ratio = w4_shots / w0_shots

    check("Shot ratio w4/w0 in [1.2, 1.6]",
          1.2 <= ratio <= 1.6,
          f"ratio={ratio:.3f}")

    check("w0 shots in [750, 900]",
          750 <= w0_shots <= 900,
          f"w0={w0_shots:.0f}")

    check("w4 shots in [1050, 1200]",
          1050 <= w4_shots <= 1200,
          f"w4={w4_shots:.0f}")

    # --- Small simulation: confirm negligible difference ---
    print("    Running quick simulation (5 trials x 5 states)...")
    n_trials = 5
    n_states_per = 5
    dep_noise = 0.03
    n_shots_uniform = 1000

    f_uniform_all = []
    f_scaled_all = []

    rng = np.random.RandomState(7777)
    for trial in range(n_trials):
        for si in range(n_states_per):
            rho = random_pure_state_rng(4, rng)
            rho_noisy = depolarizing_channel(rho, dep_noise)

            # Uniform shots
            exps_uniform = simulate_expectations(rho_noisy, ops, n_shots_uniform, rng)
            rho_u = reconstruct_robust_mle(exps_uniform, ops, n_qubits=4,
                                           n_shots=n_shots_uniform, hedge=0.05, damping=0.5)
            f_u = state_fidelity(rho_u, rho)
            f_uniform_all.append(f_u)

            # Weight-scaled shots: simulate with per-operator shot count
            exps_scaled = []
            for j, P in enumerate(ops):
                true_exp = np.trace(rho_noisy @ P).real
                ns = int(round(shots_per_op[j]))
                std = np.sqrt(max((1 - true_exp**2) / max(ns, 1), 1e-12))
                noisy_exp = true_exp + rng.normal(0, std)
                exps_scaled.append(np.clip(noisy_exp, -1, 1))

            rho_s = reconstruct_robust_mle(exps_scaled, ops, n_qubits=4,
                                           n_shots=1000, hedge=0.05, damping=0.5)
            f_s = state_fidelity(rho_s, rho)
            f_scaled_all.append(f_s)

    f_u_mean = np.mean(f_uniform_all)
    f_s_mean = np.mean(f_scaled_all)
    gap = abs(f_u_mean - f_s_mean)

    print(f"    F(uniform)={f_u_mean:.3f}, F(scaled)={f_s_mean:.3f}, |gap|={gap:.3f}")

    check("Uniform vs scaled gap < 0.05 (not significant)",
          gap < 0.05,
          f"|gap|={gap:.3f}")

    print(f"    [shot allocation took {time.time()-t0:.1f}s]")


# ====================================================================
# PART 5: MLE convergence behavior (SM Sec. XII)
# ====================================================================

def _mle_with_convergence_log(expectations, pauli_ops, n_qubits, n_shots,
                               hedge=0.05, damping=0.5, max_iter=500):
    """Run MLE and log per-iteration Frobenius norm changes."""
    dim = 2 ** n_qubits
    m = len(expectations)
    b = np.array(expectations, dtype=float)

    # Variance weighting
    if n_shots is not None and n_shots > 0:
        b_clipped = np.clip(b, -0.999, 0.999)
        variances = (1.0 - b_clipped**2) / n_shots
        weights = 1.0 / np.maximum(variances, 1e-10)
        weights /= weights.sum()
        weights *= m
    else:
        weights = np.ones(m)

    rho = np.eye(dim, dtype=complex) / dim
    P_ops = [P.copy() for P in pauli_ops]
    alpha = damping
    norms = []

    for iteration in range(max_iter):
        rho_old = rho.copy()

        predicted = np.array([np.trace(P @ rho).real for P in P_ops])

        R = hedge * np.eye(dim, dtype=complex)
        for i in range(m):
            y_i = predicted[i]
            denom = max(1.0 - y_i**2, 1e-12)
            alpha_i = (1.0 - b[i] * y_i) / denom
            beta_i = (b[i] - y_i) / denom
            R += weights[i] * (alpha_i * np.eye(dim, dtype=complex)
                               + beta_i * P_ops[i])
        R /= (m + hedge * dim)

        rho_update = R @ rho @ R.conj().T
        rho_update = (rho_update + rho_update.conj().T) / 2
        tr = np.trace(rho_update).real
        if tr > 1e-15:
            rho_update /= tr

        rho = alpha * rho_update + (1 - alpha) * rho
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > 1e-15:
            rho /= tr

        diff = np.linalg.norm(rho - rho_old, 'fro')
        norms.append(diff)

        if iteration > 20 and diff < 0.01:
            alpha = min(alpha * 1.05, 0.95)

    return rho, norms


def test_convergence():
    print("\n--- PART 5: MLE convergence behavior ---")
    t0 = time.time()

    # Load W-state hardware data (first run)
    w_file = DATA_DIR / "w_repeat_results.json"
    with open(w_file) as f:
        data = json.load(f)

    run1 = data["runs"][0]
    kstar_exps = run1["kstar_expectations"]
    rand_exps = run1["rand_expectations"]

    # Get operators
    kstar_ops, kstar_labels, kstar_indices = select_kstar_paulis(4)

    # For random arm, we need the same operators that were used
    rand_labels = run1.get("rand_labels", None)
    if rand_labels is not None:
        all_ops_full, all_labels_full, all_idx_tuples = all_pauli_operators(4)
        label_to_idx = {lbl: i for i, lbl in enumerate(all_labels_full)}
        rand_ops = [all_ops_full[label_to_idx[lbl]] for lbl in rand_labels]
    else:
        # Fallback: use random selection with same seed
        rand_ops, _, _ = select_random_paulis(4, claims.get("thm:basin", "n4_M"), seed=run1["seed"])

    # Run MLE with convergence logging
    print("    Running structured MLE (500 iters)...")
    _, norms_kstar = _mle_with_convergence_log(
        kstar_exps, kstar_ops, n_qubits=4, n_shots=1000, max_iter=500)

    print("    Running random MLE (500 iters)...")
    _, norms_rand = _mle_with_convergence_log(
        rand_exps, rand_ops, n_qubits=4, n_shots=1000, max_iter=500)

    # Analyze late-stage convergence (iterations 300-500)
    late_kstar = np.mean(norms_kstar[300:500]) if len(norms_kstar) >= 500 else np.mean(norms_kstar[-100:])
    late_rand = np.mean(norms_rand[300:500]) if len(norms_rand) >= 500 else np.mean(norms_rand[-100:])

    print(f"    Structured: mean ||drho||_F (iters 300-500) = {late_kstar:.4f}")
    print(f"    Random:     mean ||drho||_F (iters 300-500) = {late_rand:.4f}")

    # SM claims: structured ~0.03, random ~0.3
    check("Structured late ||drho||_F < 0.1",
          late_kstar < 0.1,
          f"||drho||_F={late_kstar:.4f}")

    check("Random late ||drho||_F > structured",
          late_rand > late_kstar,
          f"rand={late_rand:.4f} vs kstar={late_kstar:.4f}")

    # Check that structured converges more (ratio)
    if late_kstar > 1e-10:
        ratio = late_rand / late_kstar
        print(f"    Convergence ratio (rand/kstar) = {ratio:.1f}x")
        check("Random ||drho||_F at least 2x structured",
              ratio >= 2.0,
              f"ratio={ratio:.1f}x")
    else:
        check("Structured fully converged (||drho||_F ~ 0)",
              True, f"||drho||_F={late_kstar:.2e}")

    print(f"    [convergence took {time.time()-t0:.1f}s]")


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  CROSS-VERIFICATION: Remaining SM Claims")
    print("=" * 70)

    test_ghz_hedge_scan()      # fast, hardware data
    test_convergence()          # fast, hardware data
    test_shot_allocation()      # fast, arithmetic + small sim
    test_ablation()             # medium, simulation
    test_progressive()          # medium, simulation

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL REMAINING CLAIMS VERIFIED")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 70)
    sys.exit(0 if FAIL == 0 else 1)
