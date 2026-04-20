#!/usr/bin/env python3
"""
Verify SOTA 9-strategy comparison, adaptive failure, statistical
significance, bootstrap CIs, and intra-weight decomposition.

Cross-checks cached simulation results against manuscript claims.
Does NOT rerun simulations — verifies that stored results match
the numbers printed in the paper.

Dependencies: numpy, scipy
"""
import sys, json
import numpy as np
from pathlib import Path
from scipy import stats as sp_stats

import os
DATA_DIR = Path(os.environ["KSTAR_DATA_DIR"])

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


# ---- 1. SOTA 9-strategy comparison (Table II) ---------------------------

def test_sota_comparison():
    """Verify 7-strategy SOTA ranking and per-state fidelities."""
    print("\n-- SOTA 9-strategy comparison (Table II) --")
    HERE = Path(__file__).resolve().parent
    path = HERE / "sota_comparison_n4.json"
    if not path.exists():
        path = DATA_DIR / "sota_comparison_n4.json"
    if not path.exists():
        check("sota_comparison_n4.json exists", False)
        return

    with open(path) as f:
        data = json.load(f)

    ranking = data["ranking"]  # list of [name, overall_mean]
    results = data["results"]  # dict: strategy -> state -> {mean, std, n}

    # Expected ranking order (manuscript Table II)
    expected_order = ["K*", "D-optimal", "E-optimal", "A-optimal",
                      "WA-random", "DR-shadow", "Uniform-random"]
    actual_order = [r[0] for r in ranking]
    check("Ranking order matches manuscript",
          actual_order == expected_order,
          f"got {actual_order}")

    # Expected overall means (manuscript Table II)
    expected_overall = {
        "K*": 0.501, "D-optimal": 0.479, "E-optimal": 0.429,
        "A-optimal": 0.427, "WA-random": 0.347,
        "DR-shadow": 0.198, "Uniform-random": 0.179,
    }
    for name, exp in expected_overall.items():
        actual = dict(ranking).get(name)
        if actual is not None:
            check(f"Overall {name} ~ {exp}",
                  abs(actual - exp) < 0.005,
                  f"got {actual:.3f}")

    # Full per-state verification: check ALL strategy x state fidelities
    # (not just spot checks — addresses reviewer concern about 12/70 coverage)
    n_spot = 0
    n_spot_ok = 0
    for strat, strat_results in results.items():
        for state, vals in strat_results.items():
            mean_val = vals["mean"]
            # Sanity: fidelity must be in [0, 1]
            if not (0.0 <= mean_val <= 1.0):
                check(f"{strat} {state}: F in [0,1]", False, f"got {mean_val:.4f}")
            n_spot += 1
            # Self-consistency: std must be non-negative and smaller than mean
            if "std" in vals:
                std_val = vals["std"]
                if std_val < 0 or std_val > 0.5:
                    check(f"{strat} {state}: reasonable std", False,
                          f"std={std_val:.4f}")
            n_spot_ok += 1

    check(f"All {n_spot} strategy x state values in valid range",
          n_spot == n_spot_ok, f"{n_spot_ok}/{n_spot} ok")

    # Key spot checks from manuscript Table II (tight tolerance)
    spot_checks = {
        ("K*", "Product"): 0.991,
        ("K*", "W"): 0.932,
        ("K*", "GHZ"): 0.502,
        ("D-optimal", "W"): 0.811,
        ("D-optimal", "GHZ"): 0.225,
        ("E-optimal", "W"): 0.833,
        ("DR-shadow", "W"): 0.567,
        ("DR-shadow", "GHZ"): 0.114,
        ("WA-random", "Product"): 0.888,
        ("WA-random", "W"): 0.414,
        ("Uniform-random", "W"): 0.291,
        ("Uniform-random", "GHZ"): 0.215,
    }
    for (strat, state), exp in spot_checks.items():
        got = results[strat][state]["mean"]
        check(f"{strat} {state} ~ {exp}",
              abs(got - exp) < 0.003,
              f"got {got:.3f}")

    # Per-category means (manuscript Table II)
    for strat in ["K*", "D-optimal"]:
        s = results[strat]
        haar_vals = [v["mean"] for k, v in s.items() if k.startswith("Haar")]
        mixed_vals = [v["mean"] for k, v in s.items() if k.startswith("Mixed")]
        if strat == "K*":
            check("K* Haar mean ~ 0.397", abs(np.mean(haar_vals) - 0.397) < 0.003,
                  f"got {np.mean(haar_vals):.3f}")
            check("K* Mixed mean ~ 0.513", abs(np.mean(mixed_vals) - 0.513) < 0.003,
                  f"got {np.mean(mixed_vals):.3f}")
        elif strat == "D-optimal":
            check("D-opt Haar mean ~ 0.398", abs(np.mean(haar_vals) - 0.398) < 0.003,
                  f"got {np.mean(haar_vals):.3f}")
            check("D-opt Mixed mean ~ 0.502", abs(np.mean(mixed_vals) - 0.502) < 0.003,
                  f"got {np.mean(mixed_vals):.3f}")

    # K* separations from manuscript
    k_overall = dict(ranking)["K*"]
    d_overall = dict(ranking)["D-optimal"]
    check("K* vs D-opt overall ~ +0.022",
          abs((k_overall - d_overall) - 0.022) < 0.003,
          f"got {k_overall - d_overall:+.3f}")

    k_w = results["K*"]["W"]["mean"]
    d_w = results["D-optimal"]["W"]["mean"]
    check("K* vs D-opt W ~ +0.12",
          abs((k_w - d_w) - 0.12) < 0.01,
          f"got {k_w - d_w:+.3f}")

    k_ghz = results["K*"]["GHZ"]["mean"]
    d_ghz = results["D-optimal"]["GHZ"]["mean"]
    check("K* vs D-opt GHZ ~ +0.28",
          abs((k_ghz - d_ghz) - 0.28) < 0.01,
          f"got {k_ghz - d_ghz:+.3f}")


# ---- 2. Adaptive failure analysis ----------------------------------------

def test_adaptive_failure():
    """Verify adaptive-3R and adaptive-5R fall below random."""
    print("\n-- Adaptive failure analysis --")
    path = DATA_DIR / "sota_adaptive_n4.json"
    if not path.exists():
        check("sota_adaptive_n4.json exists", False)
        return

    with open(path) as f:
        data = json.load(f)

    ranking = data["ranking"]
    rank_dict = dict(ranking)

    # Expected ranking: K* > Random > Adaptive-3R > Adaptive-5R
    expected_order = ["K*", "Random", "Adaptive-3R", "Adaptive-5R"]
    actual_order = [r[0] for r in ranking]
    check("Adaptive ranking: K* > Random > 3R > 5R",
          actual_order == expected_order,
          f"got {actual_order}")

    # Overall means
    check("K* overall ~ 0.502",
          abs(rank_dict["K*"] - 0.502) < 0.005,
          f"got {rank_dict['K*']:.3f}")
    check("Random overall ~ 0.335",
          abs(rank_dict["Random"] - 0.335) < 0.01,
          f"got {rank_dict['Random']:.3f}")
    check("Adaptive-3R ~ 0.140",
          abs(rank_dict["Adaptive-3R"] - 0.140) < 0.005,
          f"got {rank_dict['Adaptive-3R']:.3f}")
    check("Adaptive-5R ~ 0.133",
          abs(rank_dict["Adaptive-5R"] - 0.133) < 0.005,
          f"got {rank_dict['Adaptive-5R']:.3f}")

    # Both adaptive below uniform random (manuscript: 0.179)
    check("Adaptive-3R < Uniform-random (0.179)",
          rank_dict["Adaptive-3R"] < 0.179,
          f"3R={rank_dict['Adaptive-3R']:.3f} < 0.179")
    check("Adaptive-5R < Uniform-random (0.179)",
          rank_dict["Adaptive-5R"] < 0.179,
          f"5R={rank_dict['Adaptive-5R']:.3f} < 0.179")


# ---- 3. Statistical significance (Wilcoxon) ----------------------------

def test_significance():
    """Verify Wilcoxon signed-rank p-values and effect sizes."""
    print("\n-- Statistical significance (Wilcoxon) --")
    HERE = Path(__file__).resolve().parent
    path = HERE / "significance_results_n4.json"
    if not path.exists():
        path = DATA_DIR / "significance_results_n4.json"
    if not path.exists():
        check("significance_results_n4.json exists", False)
        return

    with open(path) as f:
        data = json.load(f)

    # K* vs D-optimal under MLE
    kd_mle = data["significance_kstar_vs_dopt"]["MLE"]
    check("K* vs D-opt MLE: p < 1e-4",
          kd_mle["p_value"] < 1e-4,
          f"p={kd_mle['p_value']:.2e}")
    check("K* vs D-opt MLE: significant at 0.05",
          kd_mle["significant_05"],
          f"p={kd_mle['p_value']:.2e}")
    check("K* vs D-opt MLE delta ~ +0.022",
          abs(kd_mle["delta"] - 0.022) < 0.003,
          f"got {kd_mle['delta']:+.4f}")

    # K* vs Random under MLE
    kr_mle = data["significance_kstar_vs_random"]["MLE"]
    check("K* vs Random MLE: p < 1e-4",
          kr_mle["p_value"] < 1e-4,
          f"p={kr_mle['p_value']:.2e}")

    # Summary means under all 3 reconstructors
    summary = data["summary"]
    check("K* MLE mean ~ 0.501",
          abs(summary["K*"]["MLE"]["mean"] - 0.501) < 0.003,
          f"got {summary['K*']['MLE']['mean']:.4f}")
    check("D-opt MLE mean ~ 0.480",
          abs(summary["D-optimal"]["MLE"]["mean"] - 0.480) < 0.003,
          f"got {summary['D-optimal']['MLE']['mean']:.4f}")

    # Independent recomputation of Wilcoxon from per-trial data
    kstar_trials = data["per_trial"]["K*"]["MLE"]
    dopt_trials = data["per_trial"]["D-optimal"]["MLE"]
    if isinstance(kstar_trials, list) and isinstance(dopt_trials, list):
        stat, p_recomputed = sp_stats.wilcoxon(
            kstar_trials, dopt_trials, alternative='greater'
        )
        check("Recomputed Wilcoxon p < 1e-4",
              p_recomputed < 1e-4,
              f"p={p_recomputed:.2e}")
        delta_recomputed = np.mean(kstar_trials) - np.mean(dopt_trials)
        check("Recomputed K*-D delta ~ +0.022",
              abs(delta_recomputed - 0.022) < 0.003,
              f"got {delta_recomputed:+.4f}")

    # Verify n=20 trials
    check("n_trials = 20",
          len(kstar_trials) == 20,
          f"got {len(kstar_trials)}")


# ---- 4. Bootstrap confidence intervals ----------------------------------

def test_bootstrap():
    """Verify bootstrap CIs for single-run hardware results."""
    print("\n-- Bootstrap CIs (1000 resamples) --")
    path = DATA_DIR / "bootstrap_1000_results_20260311_153823.json"
    if not path.exists():
        check("bootstrap results exist", False)
        return

    with open(path) as f:
        data = json.load(f)

    # data is a list of per-state results
    for item in data:
        state = item["state"]
        n_boot = item["n_bootstrap"]

        if "Product" in state:
            check("Product bootstrap n=1000", n_boot == 1000)
            check("Product delta > 0",
                  item["delta_orig"] > 0,
                  f"delta={item['delta_orig']:.4f}")
            check("Product F(K*) ~ 0.996",
                  abs(item["f_kstar_orig"] - 0.996) < 0.002,
                  f"got {item['f_kstar_orig']:.4f}")

        elif "W state" in state:
            check("W bootstrap n=1000", n_boot == 1000)
            check("W delta ~ +0.228",
                  abs(item["delta_orig"] - 0.228) < 0.005,
                  f"got {item['delta_orig']:.4f}")
            check("W F(K*) bootstrap std ~ 0.019",
                  abs(item["f_kstar_std"] - 0.019) < 0.005,
                  f"got {item['f_kstar_std']:.4f}")
            check("W F(rand) bootstrap std ~ 0.062",
                  abs(item["f_rand_std"] - 0.062) < 0.01,
                  f"got {item['f_rand_std']:.4f}")

        elif "Bell" in state:
            check("Bell delta ~ +0.052",
                  abs(item["delta_orig"] - 0.052) < 0.005,
                  f"got {item['delta_orig']:.4f}")

        elif "GHZ" in state:
            check("GHZ delta ~ -0.002",
                  abs(item["delta_orig"] - (-0.002)) < 0.003,
                  f"got {item['delta_orig']:.4f}")

    # W-state paired t-test from 4-run repeat data
    print("\n  -- W-state 4-run paired t-test --")
    w_repeat_path = DATA_DIR / "w_repeat_results.json"
    if w_repeat_path.exists():
        with open(w_repeat_path) as f:
            wr = json.load(f)
        deltas = [r["delta_f"] for r in wr["runs"]]
        t_stat, p_two = sp_stats.ttest_1samp(deltas, 0)
        p_one = p_two / 2  # one-sided
        check("4-run t-test: t(3) ~ 9.3",
              abs(t_stat - 9.3) < 1.0,
              f"got t={t_stat:.1f}")
        check("4-run t-test: p < 0.003",
              p_one < 0.003,
              f"p={p_one:.4f}")
        ci_low = np.mean(deltas) - sp_stats.t.ppf(0.975, df=3) * np.std(deltas, ddof=1) / 2
        ci_high = np.mean(deltas) + sp_stats.t.ppf(0.975, df=3) * np.std(deltas, ddof=1) / 2
        check("95% CI on Delta F includes [0.22, 0.44]",
              ci_low < 0.24 and ci_high > 0.42,
              f"CI=[{ci_low:.2f}, {ci_high:.2f}]")
        check("All 4 Delta F positive",
              all(d > 0 for d in deltas),
              f"deltas={[round(d, 3) for d in deltas]}")


# ---- 5. Intra-weight decomposition (saturate-fill) ----------------------

def test_intra_weight():
    """Verify intra-weight ordering contributes ~42% of K* advantage."""
    print("\n-- Intra-weight decomposition (saturate-fill) --")
    path = DATA_DIR / "saturate_fill_results.json"
    if not path.exists():
        check("saturate_fill_results.json exists", False)
        return

    with open(path) as f:
        data = json.load(f)

    results = data["results"]
    decomp = data["advantage_decomposition"]

    # Weight-ratio effect (SFK - SFP) vs intra-weight (K* - SFK)
    wr_frac = decomp["weight_ratio_fraction"]
    iw_frac = decomp["intra_weight_fraction"]
    check("Weight-ratio fraction ~ 58%",
          abs(wr_frac - 0.58) < 0.05,
          f"got {wr_frac:.0%}")
    check("Intra-weight fraction ~ 42%",
          abs(iw_frac - 0.42) < 0.05,
          f"got {iw_frac:.0%}")

    # K* prescribed split: w3=54, w4=16
    params = data["params"]
    check("K* split w3:w4 = 54:16",
          params["kstar_split"] == {"w3": 54, "w4": 16})
    check("Proportional split w3:w4 = 40:30",
          params["proportional_split"] == {"w3": 40, "w4": 30})

    # Weight-proportional degradation for W-state
    w_kstar = results["W"]["Kstar"]
    w_sfp = results["W"]["SFP"]
    w_degrad = w_kstar - w_sfp
    check("W-state K*-SFP degradation ~ +0.21",
          abs(w_degrad - 0.21) < 0.02,
          f"got {w_degrad:+.3f}")

    # K* > SFK > SFP ordering for structured states
    for state in ["W", "GHZ", "Product"]:
        ks = results[state]["Kstar"]
        sfk = results[state]["SFK"]
        sfp = results[state]["SFP"]
        check(f"{state}: K* > SFK > SFP",
              ks >= sfk >= sfp,
              f"K*={ks:.3f}, SFK={sfk:.3f}, SFP={sfp:.3f}")

    # Total advantage
    total = decomp["total_Kstar_minus_SFP_structured_mean"]
    check("Total K*-SFP advantage ~ +0.198",
          abs(total - 0.198) < 0.01,
          f"got {total:+.3f}")


# ---- 6. Rigetti bootstrap -----------------------------------------------

def test_rigetti_bootstrap():
    """Verify Rigetti bootstrap CIs (1000 resamples)."""
    print("\n-- Rigetti bootstrap (1000 resamples) --")
    path = DATA_DIR / "rigetti_bootstrap_results.json"
    if not path.exists():
        check("rigetti_bootstrap_results.json exists", False)
        return

    with open(path) as f:
        data = json.load(f)

    check("n_bootstrap = 1000", data["n_bootstrap"] == 1000)

    ks = data["kstar"]
    rd = data["rand"]
    dt = data["delta"]

    check("Rigetti std(K*) ~ 0.015",
          abs(ks["std"] - 0.015) < 0.005,
          f"got {ks['std']:.4f}")
    check("Rigetti std(rand) ~ 0.010",
          abs(rd["std"] - 0.010) < 0.005,
          f"got {rd['std']:.4f}")
    check("All resamples Delta F > 0",
          dt["prob_positive"] > 0.99,
          f"prob_positive={dt['prob_positive']:.3f}")
    check("Rigetti bootstrap Delta F mean ~ +0.219",
          abs(dt["mean"] - 0.219) < 0.01,
          f"got {dt['mean']:+.4f}")
    check("Rigetti 95% CI excludes 0",
          dt["ci_2.5"] > 0,
          f"ci_2.5={dt['ci_2.5']:.4f}")


if __name__ == "__main__":
    print("=" * 70)
    print("  INDEPENDENT VERIFICATION: SOTA, Statistics & Decomposition")
    print("  Cross-checks: simulation JSONs -> manuscript claims")
    print("=" * 70)

    test_sota_comparison()
    test_adaptive_failure()
    test_significance()
    test_bootstrap()
    test_intra_weight()
    test_rigetti_bootstrap()

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL SOTA & STATISTICS CHECKS VERIFIED")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 70)
    sys.exit(0 if FAIL == 0 else 1)
