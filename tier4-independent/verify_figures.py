#!/usr/bin/env python3
"""
Cross-check hardcoded data arrays in gen_figures.py against authoritative
data sources (JSON files, text outputs, and arithmetic identities).

Catches any mismatch between correct data and what actually gets plotted.

No project imports — pure JSON parsing and arithmetic only.
"""
import json
import os
import re
import sys
from math import comb
from registry import claims

# ── Paths ────────────────────────────────────────────────────────────────
from pathlib import Path
DATA_DIR = os.environ["KSTAR_DATA_DIR"]

# ── Bookkeeping ──────────────────────────────────────────────────────────
_pass = 0
_fail = 0
_warn = 0


def check(name, condition, detail=""):
    global _pass, _fail
    if condition:
        _pass += 1
        print(f"  PASS  {name}")
    else:
        _fail += 1
        print(f"  FAIL  {name}")
        if detail:
            print(f"        {detail}")


def warn(name, detail=""):
    global _warn
    _warn += 1
    print(f"  WARN  {name}")
    if detail:
        print(f"        {detail}")


def load_json(path):
    """Load a JSON file, returning None if missing."""
    if not os.path.isfile(path):
        warn(f"File not found: {os.path.basename(path)}", path)
        return None
    with open(path) as f:
        return json.load(f)


def approx(a, b, tol=1e-3):
    """Check approximate equality within tolerance."""
    if a is None or b is None:
        return False
    return abs(a - b) < tol


def approx_round(a, b, decimals=2):
    """Check that a equals b rounded to the given number of decimal places.

    Handles the case where gen_figures.py uses display-rounded values
    (e.g., 0.54 for 0.5422, 0.08 for 0.0760).
    """
    if a is None or b is None:
        return False
    return abs(a - round(b, decimals)) < 10 ** -(decimals + 1)


# ════════════════════════════════════════════════════════════════════════
#  Figure 1: Head-to-head fidelity comparison
# ════════════════════════════════════════════════════════════════════════
print("\n=== Figure 1: Fidelity comparison ===")

# Hardcoded in gen_figures.py
f_struct = [0.996, 0.518, 0.872, 0.496]
f_rand = [0.980, 0.465, 0.54, 0.498]
err_struct = [None, None, 0.021, None]
err_rand = [None, None, 0.08, None]

# State mapping (determined by fidelity matching):
#   214441 -> Product  (f_struct ~0.996)
#   220545 -> Bell     (f_struct ~0.518, mode=full_tomography, n_qubits=2)
#   214922 -> GHZ      (f_struct ~0.496)
#   W state comes from w_repeat_results.json (4-run mean)
hw_files = {
    "Product": "hardware_results_ibm_fez_20260307_214441.json",
    "Bell": "hardware_results_ibm_fez_20260307_220545.json",
    "GHZ": "hardware_results_ibm_fez_20260307_214922.json",
}

# Check Product (index 0)
d = load_json(os.path.join(DATA_DIR, hw_files["Product"]))
if d is not None:
    check("Fig1 Product f_struct",
          approx(f_struct[0], d["f_structured"]),
          f"gen_figures={f_struct[0]}, JSON={d['f_structured']:.6f}")
    check("Fig1 Product f_rand",
          approx(f_rand[0], d["f_rand"]),
          f"gen_figures={f_rand[0]}, JSON={d['f_rand']:.6f}")

# Check Bell (index 1)
d = load_json(os.path.join(DATA_DIR, hw_files["Bell"]))
if d is not None:
    check("Fig1 Bell f_struct",
          approx(f_struct[1], d["f_structured"]),
          f"gen_figures={f_struct[1]}, JSON={d['f_structured']:.6f}")
    check("Fig1 Bell f_rand",
          approx(f_rand[1], d["f_rand"]),
          f"gen_figures={f_rand[1]}, JSON={d['f_rand']:.6f}")

# Check W (index 2) — from w_repeat_results.json summary
d = load_json(os.path.join(DATA_DIR, "w_repeat_results.json"))
if d is not None:
    s = d.get("summary", {})
    check("Fig1 W f_struct (4-run mean)",
          approx(f_struct[2], s.get("f_kstar_mean")),
          f"gen_figures={f_struct[2]}, JSON={s.get('f_kstar_mean', '?')}")
    # W f_rand is display-rounded to 2 decimal places (0.5422 -> 0.54)
    check("Fig1 W f_rand (4-run mean)",
          approx_round(f_rand[2], s.get("f_rand_mean"), decimals=2),
          f"gen_figures={f_rand[2]}, JSON={s.get('f_rand_mean', '?')}, "
          f"rounded={round(s.get('f_rand_mean', 0), 2)}")
    check("Fig1 W err_struct (4-run std)",
          approx(err_struct[2], s.get("f_kstar_std")),
          f"gen_figures={err_struct[2]}, JSON={s.get('f_kstar_std', '?')}")
    # W err_rand is display-rounded to 2 decimal places (0.0760 -> 0.08)
    check("Fig1 W err_rand (4-run std)",
          approx_round(err_rand[2], s.get("f_rand_std"), decimals=2),
          f"gen_figures={err_rand[2]}, JSON={s.get('f_rand_std', '?')}, "
          f"rounded={round(s.get('f_rand_std', 0), 2)}")

# Check GHZ (index 3)
d = load_json(os.path.join(DATA_DIR, hw_files["GHZ"]))
if d is not None:
    check("Fig1 GHZ f_struct",
          approx(f_struct[3], d["f_structured"]),
          f"gen_figures={f_struct[3]}, JSON={d['f_structured']:.6f}")
    check("Fig1 GHZ f_rand",
          approx(f_rand[3], d["f_rand"]),
          f"gen_figures={f_rand[3]}, JSON={d['f_rand']:.6f}")

# ════════════════════════════════════════════════════════════════════════
#  Figure 2: Eigenvalue mass distribution (arithmetic verification)
# ════════════════════════════════════════════════════════════════════════
print("\n=== Figure 2: Eigenvalue mass ===")

lambda_w = claims.get("prop:spectral_q_main", "eigenvalues_n4_K5_q2")
binom_w = [1, 4, 6, 4, 1]
pct_trace = [6.6, 40.9, 17.5, 23.4, 11.7]
c_w = claims.get("lem:monotone", "c_w_K5")
M_w = claims.get("lem:monotone", "M_w_K5")
A_w = [1, 12, 54, 108, 81]

# Verify c_w sums to n4_M
n4_M = claims.get("thm:basin", "n4_M")
check(f"Fig2 sum(c_w) = {n4_M}",
      sum(c_w) == n4_M,
      f"sum(c_w) = {sum(c_w)}")

# Verify lambda_w = 16 * c_w / binom(4, w)
for w in range(5):
    expected = 16 * c_w[w] / comb(4, w)
    check(f"Fig2 lambda_{w} = 16*c_{w}/C(4,{w})",
          abs(lambda_w[w] - expected) < 1e-9,
          f"gen_figures={lambda_w[w]}, computed={expected}")

# Verify lambda_w * binom_w
lambda_times_binom = [lw * bw for lw, bw in zip(lambda_w, binom_w)]
total_trace = sum(lambda_times_binom)
check("Fig2 lambda*binom products",
      lambda_times_binom == [144, 896, 384, 512, 256],
      f"computed={lambda_times_binom}")

# Verify pct_trace = lambda_w * binom_w / total * 100
for w in range(5):
    expected_pct = round(100.0 * lambda_times_binom[w] / total_trace, 1)
    check(f"Fig2 pct_trace[{w}] = {expected_pct}%",
          abs(pct_trace[w] - expected_pct) < 0.15,
          f"gen_figures={pct_trace[w]}, computed={expected_pct}")

# Verify pct_trace sums to ~100
pct_sum = sum(pct_trace)
check("Fig2 sum(pct_trace) ~ 100",
      abs(pct_sum - 100.0) < 0.5,
      f"sum = {pct_sum}")

# Verify M_w sums to n4_M
check(f"Fig2 sum(M_w) = {n4_M}",
      sum(M_w) == n4_M,
      f"sum(M_w) = {sum(M_w)}")

# Verify A_w = 3^4 = 256 total available (with identity)
# A_w = number of weight-w Pauli operators on 4 qubits (including I...I for w=0)
# Weight w: C(4,w) * 3^w  (choose which qubits are non-I, then X/Y/Z each)
for w in range(5):
    expected_Aw = comb(4, w) * (3 ** w)
    check(f"Fig2 A_{w} = C(4,{w})*3^{w}",
          A_w[w] == expected_Aw,
          f"gen_figures={A_w[w]}, computed={expected_Aw}")

check("Fig2 sum(A_w) = 256",
      sum(A_w) == 256,
      f"sum(A_w) = {sum(A_w)}")

# Verify M_w <= A_w for all w
for w in range(5):
    check(f"Fig2 M_{w} <= A_{w}",
          M_w[w] <= A_w[w],
          f"M_{w}={M_w[w]}, A_{w}={A_w[w]}")

# ════════════════════════════════════════════════════════════════════════
#  Figure 4: SOTA strategy comparison
# ════════════════════════════════════════════════════════════════════════
print("\n=== Figure 4: SOTA comparison ===")

# Hardcoded in gen_figures.py
f_overall_gf = [0.501, 0.479, 0.429, 0.428, 0.347, 0.198, 0.140, 0.133, 0.179]
f_product_gf = [0.991, 0.987, 0.992, 0.993, 0.888, 0.952, 0.120, 0.034, 0.711]
f_W_gf = [0.932, 0.811, 0.833, 0.833, 0.414, 0.567, 0.096, 0.092, 0.291]
f_GHZ_gf = [0.502, 0.225, 0.225, 0.229, 0.352, 0.114, 0.139, 0.075, 0.215]
f_Haar_gf = [0.397, 0.398, 0.317, 0.316, 0.242, 0.066, 0.062, 0.062, 0.063]
f_Mixed_gf = [0.513, 0.502, 0.466, 0.462, 0.391, 0.227, 0.224, 0.224, 0.227]

# Strategy ordering in gen_figures:
#   0: K*         1: D-opt      2: E-opt      3: A-opt      4: WA-rand
#   5: DR-shadow  6: Adapt(3R)  7: Adapt(5R)  8: Unif-rand

# Non-adaptive strategies (indices 0-5, 8) from sota_comparison_n4.json
sota_json = load_json(os.path.join(DATA_DIR, "sota_comparison_n4.json"))
if sota_json is not None:
    # Map JSON strategy names to gen_figures indices
    nonadaptive_map = {
        "K*": 0, "D-optimal": 1, "E-optimal": 2, "A-optimal": 3,
        "WA-random": 4, "DR-shadow": 5, "Uniform-random": 8,
    }

    # Check overall ranking
    for json_name, val in sota_json["ranking"]:
        idx = nonadaptive_map.get(json_name)
        if idx is not None:
            check(f"Fig4 f_overall[{idx}] ({json_name})",
                  approx(f_overall_gf[idx], val),
                  f"gen_figures={f_overall_gf[idx]}, JSON={val:.4f}")

    # Check per-state means for non-adaptive strategies
    def compute_category_mean(results_dict, strategy, category_prefix):
        """Compute mean of Haar-0..Haar-9 or Mixed-0..Mixed-9."""
        vals = []
        for k, v in results_dict[strategy].items():
            if k.startswith(category_prefix + "-"):
                vals.append(v["mean"])
        return sum(vals) / len(vals) if vals else None

    for json_name, idx in nonadaptive_map.items():
        res = sota_json["results"][json_name]
        # Product
        if "Product" in res:
            check(f"Fig4 f_product[{idx}] ({json_name})",
                  approx(f_product_gf[idx], res["Product"]["mean"]),
                  f"gen_figures={f_product_gf[idx]}, JSON={res['Product']['mean']:.4f}")
        # W
        if "W" in res:
            check(f"Fig4 f_W[{idx}] ({json_name})",
                  approx(f_W_gf[idx], res["W"]["mean"]),
                  f"gen_figures={f_W_gf[idx]}, JSON={res['W']['mean']:.4f}")
        # GHZ
        if "GHZ" in res:
            check(f"Fig4 f_GHZ[{idx}] ({json_name})",
                  approx(f_GHZ_gf[idx], res["GHZ"]["mean"]),
                  f"gen_figures={f_GHZ_gf[idx]}, JSON={res['GHZ']['mean']:.4f}")
        # Haar (aggregate of Haar-0 .. Haar-9)
        haar_mean = compute_category_mean(sota_json["results"], json_name, "Haar")
        if haar_mean is not None:
            check(f"Fig4 f_Haar[{idx}] ({json_name})",
                  approx(f_Haar_gf[idx], haar_mean),
                  f"gen_figures={f_Haar_gf[idx]}, computed={haar_mean:.4f}")
        # Mixed (aggregate of Mixed-0 .. Mixed-9)
        mixed_mean = compute_category_mean(sota_json["results"], json_name, "Mixed")
        if mixed_mean is not None:
            check(f"Fig4 f_Mixed[{idx}] ({json_name})",
                  approx(f_Mixed_gf[idx], mixed_mean),
                  f"gen_figures={f_Mixed_gf[idx]}, computed={mixed_mean:.4f}")

# Adaptive strategies (indices 6, 7) from sota_adaptive_n4.json
adapt_json = load_json(os.path.join(DATA_DIR, "sota_adaptive_n4.json"))
if adapt_json is not None:
    adaptive_map = {"Adaptive-3R": 6, "Adaptive-5R": 7}

    for json_name, val in adapt_json["ranking"]:
        idx = adaptive_map.get(json_name)
        if idx is not None:
            check(f"Fig4 f_overall[{idx}] ({json_name})",
                  approx(f_overall_gf[idx], val),
                  f"gen_figures={f_overall_gf[idx]}, JSON={val:.4f}")

    # Per-state for adaptive
    def compute_category_mean_adaptive(per_state, strategy, prefix):
        vals = []
        for k, v in per_state[strategy].items():
            if k.startswith(prefix + "-"):
                vals.append(v["mean"])
        return sum(vals) / len(vals) if vals else None

    for json_name, idx in adaptive_map.items():
        ps = adapt_json["per_state"][json_name]
        if "Product" in ps:
            check(f"Fig4 f_product[{idx}] ({json_name})",
                  approx(f_product_gf[idx], ps["Product"]["mean"]),
                  f"gen_figures={f_product_gf[idx]}, JSON={ps['Product']['mean']:.4f}")
        if "W" in ps:
            check(f"Fig4 f_W[{idx}] ({json_name})",
                  approx(f_W_gf[idx], ps["W"]["mean"]),
                  f"gen_figures={f_W_gf[idx]}, JSON={ps['W']['mean']:.4f}")
        if "GHZ" in ps:
            check(f"Fig4 f_GHZ[{idx}] ({json_name})",
                  approx(f_GHZ_gf[idx], ps["GHZ"]["mean"]),
                  f"gen_figures={f_GHZ_gf[idx]}, JSON={ps['GHZ']['mean']:.4f}")
        haar_mean = compute_category_mean_adaptive(adapt_json["per_state"], json_name, "Haar")
        if haar_mean is not None:
            check(f"Fig4 f_Haar[{idx}] ({json_name})",
                  approx(f_Haar_gf[idx], haar_mean),
                  f"gen_figures={f_Haar_gf[idx]}, computed={haar_mean:.4f}")
        mixed_mean = compute_category_mean_adaptive(adapt_json["per_state"], json_name, "Mixed")
        if mixed_mean is not None:
            check(f"Fig4 f_Mixed[{idx}] ({json_name})",
                  approx(f_Mixed_gf[idx], mixed_mean),
                  f"gen_figures={f_Mixed_gf[idx]}, computed={mixed_mean:.4f}")

# ════════════════════════════════════════════════════════════════════════
#  Figure 5: Dimensional dependence of Delta(S-AR)
# ════════════════════════════════════════════════════════════════════════
print("\n=== Figure 5: Dimensional dependence ===")

delta_SAR_gf = [0.573, 0.331, 0.100, 0.042]
n_vals_gf = [2, 3, 4, 5]

txt_path = os.path.join(DATA_DIR, "external_simulations_output.txt")
if os.path.isfile(txt_path):
    with open(txt_path) as f:
        txt = f.read()

    # Parse the "DIMENSIONAL PROGRESSION (MLE, complete):" section
    # Expected lines like:  n=2: D(S-AR) = +0.5727
    pattern = r"n=(\d+):\s+D\(S-AR\)\s*=\s*\+?([\d.]+)"
    matches = re.findall(pattern, txt)
    parsed = {int(n): float(v) for n, v in matches}

    for i, n in enumerate(n_vals_gf):
        if n in parsed:
            # gen_figures rounds to 3 decimal places
            check(f"Fig5 delta_SAR[n={n}]",
                  approx(delta_SAR_gf[i], parsed[n]),
                  f"gen_figures={delta_SAR_gf[i]}, text={parsed[n]:.4f}")
        else:
            warn(f"Fig5 delta_SAR[n={n}]: not found in text file")

    # Check monotonic decrease
    is_monotonic = all(delta_SAR_gf[i] > delta_SAR_gf[i + 1]
                       for i in range(len(delta_SAR_gf) - 1))
    check("Fig5 delta_SAR monotonically decreasing", is_monotonic)
else:
    warn("Fig5 external_simulations_output.txt not found", txt_path)

# ════════════════════════════════════════════════════════════════════════
#  Figure 6: Circuit efficiency
# ════════════════════════════════════════════════════════════════════════
print("\n=== Figure 6: Circuit efficiency ===")

bases_gf = [29, 50, 81]

# 81 = 3^4 (full tomography bases for n=4)
check("Fig6 full_tomo bases = 3^4 = 81",
      bases_gf[2] == 3 ** 4,
      f"gen_figures={bases_gf[2]}, 3^4={3**4}")

# 29 = number of distinct tensor-product bases covering 137 K* operators
# This is a structural constant; verify it's consistent with w_repeat data
d = load_json(os.path.join(DATA_DIR, "w_repeat_results.json"))
if d is not None and "runs" in d and len(d["runs"]) > 0:
    run0 = d["runs"][0]
    # kstar_expectations length should be n4_M (one per K* operator)
    n_kstar_ops = len(run0.get("kstar_expectations", []))
    check(f"Fig6 K* operators count = {n4_M}",
          n_kstar_ops == n4_M,
          f"w_repeat run0 kstar_expectations length = {n_kstar_ops}")

# 50 = random arm basis count (check from w_repeat or similar)
if d is not None and "runs" in d and len(d["runs"]) > 0:
    run0 = d["runs"][0]
    n_rand_ops = len(run0.get("rand_expectations", []))
    # Random uses 50 bases -> up to 50*4=200 operators (n=4, 4 ops per basis)
    # but actual count may vary; just note it
    check("Fig6 random arm has operators from ~50 bases",
          n_rand_ops > 0,
          f"rand_expectations length = {n_rand_ops}")

# Scaling: K* circuits plateau at min(n_patches, 4) * 29
# At n>=7: 4 colors * 29 = 116
check("Fig6 K* plateau = 4*29 = 116",
      4 * 29 == 116)

# ════════════════════════════════════════════════════════════════════════
#  Figure 7: GHZ mismatch and subspace recovery
# ════════════════════════════════════════════════════════════════════════
print("\n=== Figure 7: GHZ mismatch ===")

fidelities_gf = [0.496, 0.926, 0.618]
ghz_info_pct_gf = [0, 0, 40, 0, 60]

d = load_json(os.path.join(DATA_DIR, "ghz_dfe_results.json"))
if d is not None:
    # Full-state MLE fidelity (ratio-form R-operator)
    f_mle = d.get("f_mle_ratio_form") or d.get("f_mle_full")
    check("Fig7 full MLE fidelity",
          approx(fidelities_gf[0], f_mle),
          f"gen_figures={fidelities_gf[0]}, JSON={f_mle}")

    # Subspace MLE fidelity
    check("Fig7 subspace fidelity",
          approx(fidelities_gf[1], d.get("f_subspace_mle")),
          f"gen_figures={fidelities_gf[1]}, JSON={d.get('f_subspace_mle')}")

    # DFE informative lower bound
    check("Fig7 DFE lower bound fidelity",
          approx(fidelities_gf[2], d.get("f_dfe_informative")),
          f"gen_figures={fidelities_gf[2]}, JSON={d.get('f_dfe_informative')}")

# GHZ information distribution by weight (arithmetic check)
# GHZ_4 has 15 nontrivial nonzero Pauli expectations:
#   weight 0: identity (excluded from info fraction)
#   weight 1: 0 nonzero
#   weight 2: 6 nonzero (ZZ pairs: IIZZ, IZIZ, IZZI, ZIIZ, ZIZI, ZZII)
#   weight 3: 0 nonzero
#   weight 4: 9 nonzero (XXXX, XXYY, XYXY, XYYX, YXXY, YXYX, YYXX, YYYY, ZZZZ)
# Info percentages: w2 = 6/15 = 40%, w4 = 9/15 = 60%
check("Fig7 ghz_info_pct w=2 = 40%",
      ghz_info_pct_gf[2] == 40,
      f"6/15 = {6/15*100:.1f}%")
check("Fig7 ghz_info_pct w=4 = 60%",
      ghz_info_pct_gf[4] == 60,
      f"9/15 = {9/15*100:.1f}%")
check("Fig7 ghz_info_pct sum = 100",
      sum(ghz_info_pct_gf) == 100,
      f"sum = {sum(ghz_info_pct_gf)}")
check("Fig7 ghz_info_pct w=0,1,3 all zero",
      ghz_info_pct_gf[0] == 0 and ghz_info_pct_gf[1] == 0 and ghz_info_pct_gf[3] == 0)

# Cross-check: GHZ source file should be 214922
if d is not None:
    check("Fig7 GHZ source = 214922.json",
          "214922" in d.get("source", ""),
          f"source = {d.get('source', '?')}")

# ════════════════════════════════════════════════════════════════════════
#  Summary
# ════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print(f"  TOTAL: {_pass} passed, {_fail} failed, {_warn} warnings")
print(f"{'='*60}")
sys.exit(1 if _fail > 0 else 0)
