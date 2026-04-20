#!/usr/bin/env python3
"""
Verify allocation fractions and D(S-AR) values from simulation data.

Cross-checks values in external_simulations_output.txt and the
allocation percentages used in Figure 5.

Does NOT import core.py or robust_mle.py -- pure text/JSON parsing.
"""
import sys, re, json, os
import numpy as np
from pathlib import Path

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


def parse_reviewer_simulations():
    """Extract key values from external_simulations_output.txt."""
    path = DATA_DIR / "external_simulations_output.txt"
    if not path.exists():
        return None

    text = path.read_text()

    # ---- MLE data ----
    mle_data = {}

    # Format 1: tabular "n  MLE  D(S-AR)  S  AR  UR"
    pattern = r'(\d)\s+MLE\s+([\+\-]?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)'
    for m in re.finditer(pattern, text):
        n = int(m.group(1))
        mle_data[n] = {
            'D_SAR': float(m.group(2)),
            'S': float(m.group(3)),
            'AR': float(m.group(4)),
            'UR': float(m.group(5))
        }

    # Format 2: per-simulation block with "Structured: X.XXXX +/- Y.YYYY"
    sim_blocks = re.split(r'SIMULATION\s+\d+\w*:', text)
    for block in sim_blocks:
        dim_m = re.search(r'State:\s*\w+_(\d)', block)
        if not dim_m:
            dim_m = re.search(r'd=(\d)', block)
        if not dim_m:
            continue
        n = int(dim_m.group(1))
        if n in mle_data:
            continue

        s_m = re.search(r'Structured:\s+(\d+\.\d+)', block)
        ar_m = re.search(r'Alloc Random:\s+(\d+\.\d+)', block)
        ur_m = re.search(r'Uniform Random:\s+(\d+\.\d+)', block)
        dsar_m = re.search(r'Delta\(S-AR\):\s+([\+\-]?\d+\.\d+)', block)

        if s_m and ar_m and ur_m and dsar_m:
            mle_data[n] = {
                'S': float(s_m.group(1)),
                'AR': float(ar_m.group(1)),
                'UR': float(ur_m.group(1)),
                'D_SAR': float(dsar_m.group(1))
            }

    # Format 3: dimensional progression summary "n=N: D(S-AR) = +X.XXXX"
    for m in re.finditer(r'n=(\d):\s*D\(S-AR\)\s*=\s*([\+\-]?\d+\.\d+)', text):
        n = int(m.group(1))
        if n not in mle_data:
            mle_data[n] = {'D_SAR': float(m.group(2))}
        elif 'D_SAR' not in mle_data[n]:
            mle_data[n]['D_SAR'] = float(m.group(2))

    # ---- PLS data (tabular format) ----
    pls_data = {}
    pattern_pls = r'(\d)\s+PLS\s+([\+\-]?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)'
    for m in re.finditer(pattern_pls, text):
        n = int(m.group(1))
        pls_data[n] = {'D_SAR': float(m.group(2))}

    return mle_data, pls_data


def test_dsar_values():
    """Verify D(S-AR) values match manuscript."""
    print("\n-- D(S-AR) progression (MLE) --")
    result = parse_reviewer_simulations()
    if result is None:
        print("  [SKIP] external_simulations_output.txt not in data archive (simulation-only)")
        return

    mle_data, pls_data = result

    expected_dsar = {2: 0.573, 3: 0.331, 4: 0.100, 5: 0.042}
    for n, exp in expected_dsar.items():
        if n in mle_data and 'D_SAR' in mle_data[n]:
            got = mle_data[n]['D_SAR']
            check(f"n={n}: D(S-AR) ~ {exp}",
                  abs(got - exp) < 0.002,
                  f"data={got:.4f}, manuscript={exp}")
        else:
            check(f"n={n}: data found", False)


def test_allocation_fractions():
    """Verify allocation fraction = D(AR-UR) / D(S-UR)."""
    print("\n-- Allocation fractions --")
    result = parse_reviewer_simulations()
    if result is None:
        print("  [SKIP] external_simulations_output.txt not in data archive (simulation-only)")
        return

    mle_data, _ = result

    # n=2 through n=5 (n=5 allocation ~86%, "approaches unity")
    expected_alloc = {2: -13, 3: 43, 4: 81, 5: 86}
    for n, exp_pct in expected_alloc.items():
        if n not in mle_data or 'S' not in mle_data[n]:
            check(f"n={n}: full S/AR/UR data found", False)
            continue

        d = mle_data[n]
        d_ar_ur = d['AR'] - d['UR']
        d_s_ur = d['S'] - d['UR']

        if abs(d_s_ur) < 1e-6:
            check(f"n={n}: D(S-UR) nonzero", False)
            continue

        alloc_frac = 100 * d_ar_ur / d_s_ur
        check(f"n={n}: alloc ~ {exp_pct}%",
              abs(alloc_frac - exp_pct) < 3,
              f"D(AR-UR)={d_ar_ur:.4f}, D(S-UR)={d_s_ur:.4f}, frac={alloc_frac:.1f}%")

    # Specific Delta(AR-UR) value at n=4 (manuscript: +0.431)
    if 4 in mle_data and 'AR' in mle_data[4] and 'UR' in mle_data[4]:
        d_ar_ur_4 = mle_data[4]['AR'] - mle_data[4]['UR']
        check("n=4: D(AR-UR) ~ +0.431",
              abs(d_ar_ur_4 - 0.431) < 0.005,
              f"got {d_ar_ur_4:.4f}")


def test_mle_ratio():
    """Verify 13.8x MLE ratio and 28.3x PLS ratio."""
    print("\n-- Dimensional ratios --")
    result = parse_reviewer_simulations()
    if result is None:
        return

    mle_data, pls_data = result

    if 2 in mle_data and 5 in mle_data:
        d2 = mle_data[2]['D_SAR']
        d5 = mle_data[5]['D_SAR']
        mle_ratio = d2 / d5
        check("MLE ratio n=2/n=5 ~ 13.8",
              abs(mle_ratio - 13.8) < 0.2,
              f"{d2:.4f}/{d5:.4f} = {mle_ratio:.2f}")

    if 2 in pls_data and 5 in pls_data:
        p2 = pls_data[2]['D_SAR']
        p5 = pls_data[5]['D_SAR']
        pls_ratio = p2 / p5
        check("PLS ratio n=2/n=5 ~ 28.3",
              abs(pls_ratio - 28.3) < 0.2,
              f"{p2}/{p5} = {pls_ratio:.2f}")


def test_three_arm_hardware():
    """Verify three-arm hardware results and allocation fraction."""
    print("\n-- Three-arm hardware (ibm_fez) --")
    path = DATA_DIR / "three_arm_hardware_ibm_fez_20260311_142656.json"
    if not path.exists():
        check("three-arm data exists", False)
        return

    with open(path) as f:
        data = json.load(f)

    runs = data.get("runs", data.get("results", []))
    f_s = [r["f_kstar"] for r in runs]
    f_ar = [r["f_ar"] for r in runs]
    f_ur = [r["f_ur"] for r in runs]

    mean_s = np.mean(f_s)
    mean_ar = np.mean(f_ar)
    mean_ur = np.mean(f_ur)

    print(f"  F(S) = {mean_s:.3f}, F(AR) = {mean_ar:.3f}, F(UR) = {mean_ur:.3f}")

    # All three fidelity means
    check("F(S) ~ 0.903", abs(mean_s - 0.903) < 0.02, f"got {mean_s:.3f}")
    check("F(AR) ~ 0.722", abs(mean_ar - 0.722) < 0.02, f"got {mean_ar:.3f}")
    check("F(UR) ~ 0.256", abs(mean_ur - 0.256) < 0.05, f"got {mean_ur:.3f}")

    # Allocation fraction
    d_ar_ur = mean_ar - mean_ur
    d_s_ur = mean_s - mean_ur
    if abs(d_s_ur) > 1e-6:
        alloc = 100 * d_ar_ur / d_s_ur
        check("Hardware alloc ~ 76%", abs(alloc - 76) < 15,
              f"got {alloc:.0f}%")

    # Hardware Delta(S-AR) = +0.181 mean
    d_s_ar = mean_s - mean_ar
    check("Hardware D(S-AR) ~ +0.181", abs(d_s_ar - 0.181) < 0.02,
          f"got {d_s_ar:.3f}")

    # Per-seed values (spot check)
    n_seeds = len(runs)
    check(f"Three seeds measured", n_seeds == 3, f"got {n_seeds}")

    # All seeds: F(S) > F(AR) > F(UR) ordering
    ordering_ok = all(r["f_kstar"] > r["f_ar"] for r in runs)
    check("All seeds: F(K*) > F(AR)", ordering_ok)


if __name__ == "__main__":
    print("=" * 70)
    print("  INDEPENDENT VERIFICATION: Allocation & D(S-AR)")
    print("  Cross-checks: simulation data -> manuscript claims")
    print("=" * 70)

    test_dsar_values()
    test_allocation_fractions()
    test_mle_ratio()
    test_three_arm_hardware()

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL ALLOCATION CHECKS VERIFIED")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 70)
    sys.exit(0 if FAIL == 0 else 1)
