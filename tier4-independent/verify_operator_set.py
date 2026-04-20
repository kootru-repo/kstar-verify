#!/usr/bin/env python3
"""
Verify that the 137 K*-selected Pauli operators match the lattice theory.

Independently enumerates lattice vectors in Z^4 with |m|^2 <= 5,
maps them to Pauli weight classes via parity reduction m -> m mod 2,
and checks that the project's select_kstar_paulis(4) returns operators
consistent with this weight budget.

Uses: numpy, math, sys.path to import core.py (project code under test).
"""
import sys, json, math
import numpy as np
from itertools import product as cart_product
from pathlib import Path

import os
# Import the project's core.py
from core import select_kstar_paulis, select_random_paulis
from registry import claims

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


def pauli_weight(label):
    """Hamming weight of Pauli label (count of non-I characters)."""
    return sum(1 for c in label if c != 'I')


def greedy_basis_cover(labels, n_qubits):
    """Greedy set cover: find minimal tensor-product bases for Pauli labels.

    A tensor-product basis is a string in {X,Y,Z}^n.
    Operator P is compatible with basis B if every non-I position matches.
    Returns list of basis strings (deterministic: sorted + lexicographic tie-break).
    """
    def compatible_bases(label):
        choices = []
        for c in label:
            choices.append(['X', 'Y', 'Z'] if c == 'I' else [c])
        return [''.join(b) for b in cart_product(*choices)]

    def covers(basis, label):
        return all(label[i] == 'I' or label[i] == basis[i]
                   for i in range(n_qubits))

    uncovered = set(labels) - {'I' * n_qubits}
    bases = []
    while uncovered:
        best_basis, best_cover = None, set()
        candidates = set()
        for lbl in uncovered:
            for b in compatible_bases(lbl):
                candidates.add(b)
        for b in sorted(candidates):
            covered = {lbl for lbl in uncovered if covers(b, lbl)}
            if len(covered) > len(best_cover) or (
                len(covered) == len(best_cover) and
                (best_basis is None or b < best_basis)
            ):
                best_basis, best_cover = b, covered
        bases.append(best_basis)
        uncovered -= best_cover
    return bases


def test_operator_count_and_weights():
    """K* selection returns 137 operators with correct weight distribution."""
    print("\n-- K* operator set (n=4) --")
    ops, labels, indices = select_kstar_paulis(4)

    n4_M = claims.get("thm:basin", "n4_M")
    check(f"len(operators) = {n4_M}", len(ops) == n4_M, f"got {len(ops)}")
    check(f"len(labels) = {n4_M}", len(labels) == n4_M)
    check(f"len(indices) = {n4_M}", len(indices) == n4_M)

    # Weight distribution
    weight_counts = [0] * 5
    for lbl in labels:
        weight_counts[pauli_weight(lbl)] += 1

    expected_M = claims.get("lem:monotone", "M_w_K5")
    for w in range(5):
        check(f"M_{w} = {expected_M[w]}", weight_counts[w] == expected_M[w],
              f"got {weight_counts[w]}")

    # Identity must be included
    check("IIII in K* set", "IIII" in labels)

    # Labels are valid Pauli strings
    valid = all(len(l) == 4 and all(c in "IXYZ" for c in l) for l in labels)
    check("All labels are valid 4-qubit Pauli strings", valid)

    # No duplicates
    check("No duplicate labels", len(set(labels)) == len(labels))

    # Saturation check: weights 0-2 use ALL available Paulis
    A_w = [math.comb(4, w) * 3**w for w in range(5)]
    for w in range(3):
        check(f"Weight {w} saturated (M_{w} = A_{w})",
              weight_counts[w] == A_w[w],
              f"M={weight_counts[w]}, A={A_w[w]}")


def test_weight_budget_from_lattice():
    """Independently derive weight budget from Z^4 lattice — no core.py.

    The allocation rule is: M_w = min(c_w, A_w), where c_w counts lattice
    points by parity weight and A_w = C(n,w)*3^w is the available Pauli pool.
    This is the Krawtchouk eigenvalue-mass allocation, derived from first
    principles without importing any project code.
    """
    print("\n-- Lattice-derived weight budget (standalone derivation) --")
    d, K = 4, 5

    # Step 1: Count lattice vectors by parity weight (brute force Z^4)
    bound = math.isqrt(K) + 1
    c_w = [0] * (d + 1)
    for m in cart_product(range(-bound, bound + 1), repeat=d):
        if sum(x*x for x in m) <= K:
            hw = sum(1 for x in m if x % 2 != 0)
            c_w[hw] += 1

    expected_c = claims.get("lem:monotone", "c_w_K5")
    for w in range(5):
        check(f"c_{w} = {expected_c[w]}", c_w[w] == expected_c[w], f"got {c_w[w]}")

    n4_M = claims.get("thm:basin", "n4_M")
    check(f"sum(c_w) = {n4_M}", sum(c_w) == n4_M, f"got {sum(c_w)}")

    # Step 2: Available Paulis per weight
    A_w = [math.comb(4, w) * 3**w for w in range(5)]
    expected_A = [1, 12, 54, 108, 81]
    for w in range(5):
        check(f"A_{w} = C(4,{w})*3^{w} = {expected_A[w]}", A_w[w] == expected_A[w])

    # Step 3: DERIVE M_w via bottom-up saturation (no hardcoded values)
    #   (a) Initial: M_w = min(c_w, A_w) — lattice budget capped at available
    #   (b) Collect surplus from weights where c_w > A_w
    #   (c) Redistribute surplus bottom-up to fill unsaturated weights
    # This follows from the Krawtchouk eigenvalue-mass principle: low-weight
    # operators carry the most information per operator, so saturate first.
    M_w_derived = [min(c_w[w], A_w[w]) for w in range(5)]
    surplus = sum(max(0, c_w[w] - A_w[w]) for w in range(5))
    check(f"Initial allocation = {sum(M_w_derived)}, surplus = {surplus}",
          sum(M_w_derived) + surplus == sum(c_w),
          f"{sum(M_w_derived)} + {surplus} = {sum(M_w_derived) + surplus}")

    for w in range(5):
        gap = A_w[w] - M_w_derived[w]
        if gap > 0 and surplus > 0:
            fill = min(gap, surplus)
            M_w_derived[w] += fill
            surplus -= fill

    expected_M = claims.get("lem:monotone", "M_w_K5")
    for w in range(5):
        check(f"Derived M_{w} = {expected_M[w]}",
              M_w_derived[w] == expected_M[w],
              f"got {M_w_derived[w]}")

    n4_M = claims.get("thm:basin", "n4_M")
    check(f"sum(M_w) = {n4_M}", sum(M_w_derived) == n4_M, f"got {sum(M_w_derived)}")
    check("surplus fully distributed", surplus == 0, f"remaining surplus = {surplus}")

    # Step 4: Verify core.py's selection matches the derived allocation
    _, labels, _ = select_kstar_paulis(4)
    weight_counts = [0] * 5
    for lbl in labels:
        weight_counts[pauli_weight(lbl)] += 1
    for w in range(5):
        check(f"core.py M_{w} matches derived",
              weight_counts[w] == M_w_derived[w],
              f"core={weight_counts[w]}, derived={M_w_derived[w]}")

    # Step 5: Verify saturation — weights 0-2 use ALL available Paulis
    w_sat = -1
    for w in range(5):
        if M_w_derived[w] == A_w[w]:
            w_sat = w
        else:
            break
    check("w_sat = 2 (derived)", w_sat == 2, f"got {w_sat}")


def test_basis_compression():
    """Verify 137 K* operators compress into 29 tensor-product bases."""
    print("\n-- Basis compression --")
    _, labels, _ = select_kstar_paulis(4)

    bases = greedy_basis_cover(labels, 4)
    n_bases = len(bases)
    check("K* greedy cover = 29 bases", n_bases == 29, f"got {n_bases}")

    # K* operators touch all 81 tensor-product bases (paper claim)
    all_81 = set(''.join(b) for b in cart_product('XYZ', repeat=4))
    touched = set()
    for lbl in labels:
        for b in all_81:
            if all(lbl[i] == 'I' or lbl[i] == b[i] for i in range(4)):
                touched.add(b)
    check("K* operators touch all 81 bases", len(touched) == 81,
          f"touched {len(touched)}/81")

    # Random selection compression
    _, rand_labels, _ = select_random_paulis(4, claims.get("thm:basin", "n4_M"), seed=42)
    bases_r = greedy_basis_cover(rand_labels, 4)
    n_bases_r = len(bases_r)
    check("Random compresses to ~50 bases", 45 <= n_bases_r <= 55, f"got {n_bases_r}")
    ratio = round(n_bases_r / n_bases, 1)
    check("Compression ratio ~ 1.7x", abs(ratio - 1.7) < 0.2, f"got {ratio}x")


def test_rigetti_operator_match():
    """Verify Rigetti data has same operators as K* selection."""
    print("\n-- Rigetti operator match --")
    rig_file = DATA_DIR / "oq_grouped_results_20260316.json"
    if not rig_file.exists():
        check("Rigetti data file exists", False, f"not found: {rig_file}")
        return

    with open(rig_file) as f:
        data = json.load(f)

    kstar_exp = data["results"]["n4_kstar"]["expectations"]
    n_ops_rig = len(kstar_exp)
    n4_M = claims.get("thm:basin", "n4_M")
    check(f"Rigetti K* has {n4_M} operators", n_ops_rig == n4_M, f"got {n_ops_rig}")

    _, labels, _ = select_kstar_paulis(4)
    kstar_labels_set = set(labels)
    rig_labels_set = set(kstar_exp.keys())
    check("Rigetti K* labels match core.py", kstar_labels_set == rig_labels_set,
          f"missing from Rigetti: {kstar_labels_set - rig_labels_set}, "
          f"extra in Rigetti: {rig_labels_set - kstar_labels_set}")

    check("Rigetti n_bases = 29", data["results"]["n4_kstar"]["n_bases"] == 29)
    check(f"Rigetti n_operators = {n4_M}", data["results"]["n4_kstar"]["n_operators"] == n4_M)


if __name__ == "__main__":
    print("=" * 70)
    print("  INDEPENDENT VERIFICATION: K* Operator Set")
    print("  Checks: weight budget, lattice correspondence,")
    print("  basis compression, Rigetti operator match")
    print("=" * 70)

    test_operator_count_and_weights()
    test_weight_budget_from_lattice()
    test_basis_compression()
    test_rigetti_operator_match()

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed")
    if FAIL == 0:
        print("  ALL OPERATOR SET CHECKS VERIFIED")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 70)
    sys.exit(0 if FAIL == 0 else 1)
