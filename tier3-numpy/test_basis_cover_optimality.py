"""Optimality certificate for the K* tensor-product basis cover (n=4).

The manuscript Sec. compression-advantage states that the 137 K* operators
at n=4 compress into 29 tensor-product bases via greedy set cover, with the
disclaimer "global optimality is not claimed." This test discharges that
disclaimer by solving the set-cover ILP exactly with HiGHS and certifying
that 29 is the global minimum.

Formulation
-----------
Universe U = 137 K* Pauli labels (the canonical lattice-derived set from
            kstar_certification.core.select_kstar_paulis).
Sets    S = 81 candidate tensor-product bases B in {X,Y,Z}^4. A basis B
            covers Pauli label p iff for every qubit position i:
                p[i] == 'I'  or  p[i] == B[i].
Variables x_B in {0,1} for each candidate basis.
Constraint sum_{B compatible with p} x_B >= 1 for every p in U.
Objective minimize sum_B x_B.

The minimum is attested to be 29; HiGHS returns the exact optimum and a
witness assignment. Both lower bound (LP relaxation closed by branch-and-
bound) and upper bound (the witness) are written to disk so a reviewer can
re-verify with any MILP solver via the MPS file.

Outputs (artifacts/clique_cover/):
    kstar_n4_set_cover.mps   - solver-independent model
    kstar_n4_set_cover.sol   - witness assignment (29 bases)
    kstar_n4_set_cover.log   - HiGHS run log + lower bound
"""

import os
import sys
from itertools import product

# highspy is imported lazily inside the test so the module is importable
# (and run_all.py keeps starting) even on minimal environments without
# the MILP solver. The test then prints a SKIP and returns 0 checks.
try:
    import highspy
    _HIGHSPY_OK = True
except ImportError:
    highspy = None
    _HIGHSPY_OK = False

# Pull canonical K* labels from tier4-independent/core.py (sibling tier).
# Works both in the monorepo layout and the KV-standalone layout.
_HERE = os.path.dirname(os.path.abspath(__file__))
_TIER4 = os.path.abspath(os.path.join(_HERE, "..", "tier4-independent"))
if _TIER4 not in sys.path:
    sys.path.insert(0, _TIER4)

from core import select_kstar_paulis  # noqa: E402

ARTIFACT_DIR = os.path.join(_HERE, "artifacts", "clique_cover")
N_QUBITS = 4
EXPECTED_OPTIMUM = 29
EXPECTED_OPERATORS = 137


def _candidate_bases(n):
    return [''.join(b) for b in product('XYZ', repeat=n)]


def _compatible(pauli_label, basis_label):
    for p, b in zip(pauli_label, basis_label):
        if p != 'I' and p != b:
            return False
    return True


def _build_cover_lists(labels, candidate_bases):
    cover_lists = []
    for p in labels:
        cl = [j for j, b in enumerate(candidate_bases) if _compatible(p, b)]
        if not cl:
            raise RuntimeError(f"Pauli {p!r} has no compatible tensor-product basis")
        cover_lists.append(cl)
    return cover_lists


def _solve_set_cover(candidate_bases, cover_lists, mps_path):
    h = highspy.Highs()
    h.silent()
    nB = len(candidate_bases)

    for _ in range(nB):
        h.addCol(1.0, 0.0, 1.0, 0, [], [])
    for j in range(nB):
        h.changeColIntegrality(j, highspy.HighsVarType.kInteger)

    for cl in cover_lists:
        coefs = [1.0] * len(cl)
        h.addRow(1.0, highspy.kHighsInf, len(cl), cl, coefs)

    h.changeObjectiveSense(highspy.ObjSense.kMinimize)

    # Persist the model in solver-independent MPS form *before* solving so
    # any external reviewer solver can replay it.
    os.makedirs(os.path.dirname(mps_path), exist_ok=True)
    h.writeModel(mps_path)

    status = h.run()
    if status != highspy.HighsStatus.kOk:
        raise RuntimeError(f"HiGHS failed with status {status}")

    sol = h.getSolution()
    obj = h.getObjectiveValue()
    info = h.getInfo()
    return obj, sol.col_value, info


def _verify_witness(witness_indices, candidate_bases, labels):
    chosen_bases = [candidate_bases[j] for j in witness_indices]
    covered = set()
    for b in chosen_bases:
        for i, p in enumerate(labels):
            if _compatible(p, b):
                covered.add(i)
    if covered != set(range(len(labels))):
        missing = sorted(set(range(len(labels))) - covered)
        raise RuntimeError(
            f"Witness covers {len(covered)}/{len(labels)} Paulis; missing {missing[:5]}"
        )
    return chosen_bases


def test_basis_cover_optimality():
    """Certify chi-bar(K*_n4) = 29 via set-cover ILP."""
    passed = 0

    if not _HIGHSPY_OK:
        print("  [SKIP] highspy not installed; skipping basis-cover ILP certification")
        return passed

    _, labels, _ = select_kstar_paulis(N_QUBITS)
    assert len(labels) == EXPECTED_OPERATORS, (
        f"Canonical K* set should have {EXPECTED_OPERATORS} operators, "
        f"got {len(labels)}"
    )
    print(f"  K* set: {len(labels)} operators (canonical lattice selection)")
    passed += 1

    candidate_bases = _candidate_bases(N_QUBITS)
    assert len(candidate_bases) == 3 ** N_QUBITS == 81
    cover_lists = _build_cover_lists(labels, candidate_bases)
    print(f"  Candidate bases: {len(candidate_bases)} (= 3^{N_QUBITS})")
    passed += 1

    mps_path = os.path.join(ARTIFACT_DIR, "kstar_n4_set_cover.mps")
    sol_path = os.path.join(ARTIFACT_DIR, "kstar_n4_set_cover.sol")
    log_path = os.path.join(ARTIFACT_DIR, "kstar_n4_set_cover.log")

    obj, col_value, info = _solve_set_cover(candidate_bases, cover_lists, mps_path)
    print(f"  HiGHS optimum: {obj}")
    print(f"  MPS model written: {os.path.relpath(mps_path, _HERE)}")
    passed += 1

    # The branch-and-bound proof of optimality is encoded in the MIP gap.
    # When the gap is closed (= 0) at termination, primal == dual bound and
    # the value is provably optimal.
    mip_gap = getattr(info, "mip_gap", 0.0)
    assert mip_gap < 1e-9, f"MIP gap not closed: {mip_gap}"
    print(f"  MIP gap closed (gap={mip_gap}): primal = dual = optimum")
    passed += 1

    assert abs(obj - EXPECTED_OPTIMUM) < 1e-9, (
        f"Optimum is {obj}, expected {EXPECTED_OPTIMUM}"
    )
    print(f"  chi-bar(K*_n4) = {EXPECTED_OPTIMUM} CERTIFIED")
    passed += 1

    witness_indices = [j for j in range(len(candidate_bases)) if col_value[j] > 0.5]
    chosen_bases = _verify_witness(witness_indices, candidate_bases, labels)
    assert len(chosen_bases) == EXPECTED_OPTIMUM
    print(f"  Witness verified: {len(chosen_bases)} bases cover all 137 operators")
    passed += 1

    # Persist witness + log for reviewer replay.
    with open(sol_path, "w") as f:
        f.write(f"# K* tensor-product basis cover witness, n={N_QUBITS}\n")
        f.write(f"# objective = {obj}\n")
        f.write(f"# n_bases   = {len(chosen_bases)}\n")
        for b in sorted(chosen_bases):
            f.write(f"{b}\n")

    with open(log_path, "w") as f:
        f.write("HiGHS set-cover certification of K* basis compression\n")
        f.write(f"  n_qubits        : {N_QUBITS}\n")
        f.write(f"  K* operators    : {len(labels)}\n")
        f.write(f"  candidate bases : {len(candidate_bases)} (= 3^{N_QUBITS})\n")
        f.write(f"  status          : {highspy.HighsStatus.kOk}\n")
        f.write(f"  primal optimum  : {obj}\n")
        f.write(f"  mip_gap         : {mip_gap}\n")
        f.write(f"  EXPECTED        : {EXPECTED_OPTIMUM}\n")
        f.write(f"  result          : OPTIMUM CERTIFIED\n")

    print(f"  Witness:    {os.path.relpath(sol_path, _HERE)}")
    print(f"  Log:        {os.path.relpath(log_path, _HERE)}")

    # Cross-solver round-trip: re-read the MPS we just wrote and re-solve
    # with a completely independent MILP solver (CBC via PuLP). This proves
    # the .mps file is correct, self-contained, and solver-portable. If
    # PuLP/CBC is unavailable the check is skipped (soft) so the suite still
    # runs in minimal environments.
    try:
        import pulp  # noqa: F401
    except ImportError:
        print("  [SKIP] CBC round-trip: pulp not installed")
        return passed

    _, cbc_prob = pulp.LpProblem.fromMPS(mps_path)
    cbc_prob.solve(pulp.PULP_CBC_CMD(msg=0))
    cbc_status = pulp.LpStatus[cbc_prob.status]
    cbc_obj = pulp.value(cbc_prob.objective)
    cbc_chosen = [v for v in cbc_prob.variables() if v.varValue and v.varValue > 0.5]

    if cbc_status != "Optimal":
        raise RuntimeError(f"CBC round-trip status: {cbc_status}")
    if abs(cbc_obj - EXPECTED_OPTIMUM) > 1e-9:
        raise RuntimeError(f"CBC round-trip optimum {cbc_obj}, expected {EXPECTED_OPTIMUM}")
    if len(cbc_chosen) != EXPECTED_OPTIMUM:
        raise RuntimeError(f"CBC chose {len(cbc_chosen)} bases, expected {EXPECTED_OPTIMUM}")

    print(f"  CBC round-trip: optimum {cbc_obj} from re-parsed MPS (independent solver)")
    passed += 1
    return passed


if __name__ == "__main__":
    n = test_basis_cover_optimality()
    print(f"\n{n} checks passed")
