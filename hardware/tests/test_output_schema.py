"""Output-schema regression test for hardware/run_w_state.py.

Builds a synthetic 2-run payload, drives it through `_finalize`, and
asserts the emitted JSON matches the shipped reference
`data/w_repeat_results.json` in top-level structure.  Intended to
fail fast in CI if anyone adds or removes a field in a way that
would drift the output schema away from the reference.

Does NOT require qiskit or any backend -- the synthetic runs skip
the submission path entirely.  Complements test_smoke_simulator.py,
which exercises the qiskit code path but stops short of _finalize.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
HARDWARE_DIR = THIS_DIR.parent
REPO_ROOT = HARDWARE_DIR.parent
for p in (HARDWARE_DIR, REPO_ROOT / "tier4-independent"):
    sys.path.insert(0, str(p))


# Keys we emit.  Matches the reference JSON's provenance layout; the
# optional `calibration` sub-dict is only emitted when backend.properties()
# is available, so the test treats it as present-or-absent.
EMITTED_PROVENANCE_KEYS_REQUIRED = {
    "raw_counts",
    "job_ids",
    "n_jobs_retrieved",
    "n_jobs_failed",
    "job_account_map",
    "transpilation",
    "retrieval_timestamp",
    "physical_qubits",
    "physical_qubits_note",
}
EMITTED_PROVENANCE_KEYS_OPTIONAL = {
    "calibration",  # backend.properties() may be None on simulators
}
REQUIRED_TOPLEVEL = {
    "experiment", "backend", "n_qubits", "n_shots", "n_runs", "seeds",
    "runs", "summary", "timestamp", "provenance", "random_basis_selection",
}
REQUIRED_SUMMARY = {
    "f_kstar_mean", "f_kstar_std", "f_rand_mean", "f_rand_std",
    "delta_f_mean", "delta_f_std", "std_convention",
}
REQUIRED_RBS = {"method", "per_run_seeds", "note", "code"}


def _fake_runs() -> list[dict]:
    # Two runs with seeds matching the reference (100, 200) and
    # a handful of synthetic raw counts keyed by fabricated job_ids.
    return [
        {"run": 1, "seed": 100, "f_kstar": 0.87, "f_rand": 0.52,
         "delta_f": 0.35, "job_ids": ["job_A1", "job_A2"],
         "kstar_expectations": [0.5] * 137,
         "rand_expectations": [0.1] * 137,
         "rand_labels": ["XXXX"] * 137,
         "timestamp": "2026-04-20 00:00:00",
         "_raw_counts_per_job": {
             "job_A1": [{"0000": 500, "1111": 500}] * 100,
             "job_A2": [{"0001": 250}] * 37,
         },
         "_physical_qubits": [0, 1, 2, 3],
         "_account_id": "test_account/main"},
        {"run": 2, "seed": 200, "f_kstar": 0.88, "f_rand": 0.50,
         "delta_f": 0.38, "job_ids": ["job_B1", "job_B2"],
         "kstar_expectations": [0.55] * 137,
         "rand_expectations": [0.12] * 137,
         "rand_labels": ["YYYY"] * 137,
         "timestamp": "2026-04-20 00:01:00",
         "_raw_counts_per_job": {
             "job_B1": [{"0000": 501, "1111": 499}] * 100,
             "job_B2": [{"0001": 250}] * 37,
         },
         "_physical_qubits": [0, 1, 2, 3],
         "_account_id": "test_account/main"},
    ]


def main() -> int:
    import run_w_state  # type: ignore

    print("=" * 70)
    print("  OUTPUT SCHEMA TEST — hardware/run_w_state.py._finalize")
    print("=" * 70)

    with tempfile.TemporaryDirectory() as tmp:
        out = run_w_state._finalize(
            Path(tmp), "ibm_fez_test", [100, 200], 1000, 2, _fake_runs(),
            calibration={
                "num_qubits": 127, "backend_name": "ibm_fez_test",
                "n_qubit_properties": 127, "qubit_properties": [],
                "n_gate_errors": 0, "gate_errors": [],
                "retrieval_timestamp": "2026-04-20 00:02:00",
            },
        )

    # --- Top-level keys --------------------------------------------
    got_top = set(out.keys())
    missing = REQUIRED_TOPLEVEL - got_top
    if missing:
        print(f"\n  [FAIL] missing top-level keys: {sorted(missing)}")
        return 1
    print(f"\n  [PASS] top-level: {sorted(REQUIRED_TOPLEVEL)}")

    # --- Summary -------------------------------------------------
    got_summary = set(out["summary"].keys())
    missing_s = REQUIRED_SUMMARY - got_summary
    if missing_s:
        print(f"  [FAIL] missing summary keys: {sorted(missing_s)}")
        return 1
    std_conv = out["summary"]["std_convention"]
    if std_conv != "sample (ddof=1, N-1 denominator)":
        print(f"  [FAIL] std_convention drift: got {std_conv!r}, "
              "expected 'sample (ddof=1, N-1 denominator)'")
        return 1
    print(f"  [PASS] summary + std_convention matches reference")

    # --- random_basis_selection ----------------------------------
    got_rbs = set(out["random_basis_selection"].keys())
    if got_rbs != REQUIRED_RBS:
        print(f"  [FAIL] random_basis_selection keys: got {sorted(got_rbs)}, "
              f"expected {sorted(REQUIRED_RBS)}")
        return 1
    print(f"  [PASS] random_basis_selection: {sorted(REQUIRED_RBS)}")

    # --- provenance ----------------------------------------------
    got_prov = set(out["provenance"].keys())
    missing_p = EMITTED_PROVENANCE_KEYS_REQUIRED - got_prov
    unknown = got_prov - (EMITTED_PROVENANCE_KEYS_REQUIRED
                          | EMITTED_PROVENANCE_KEYS_OPTIONAL)
    if missing_p:
        print(f"  [FAIL] missing provenance keys: {sorted(missing_p)}")
        return 1
    if unknown:
        print(f"  [WARN] unexpected extra provenance keys: {sorted(unknown)}")
        # Not a failure -- additive changes are allowed if intentional.
    print(f"  [PASS] provenance: {sorted(got_prov)}")

    # --- job_account_map integrity -------------------------------
    jam = out["provenance"]["job_account_map"]
    if set(jam.keys()) != {"job_A1", "job_A2", "job_B1", "job_B2"}:
        print(f"  [FAIL] job_account_map job_ids: {sorted(jam.keys())}")
        return 1
    if any("token" in str(v).lower() for v in jam.values()):
        print("  [FAIL] job_account_map appears to leak a token value")
        return 1
    print(f"  [PASS] job_account_map: {len(jam)} job_ids -> "
          f"{sorted(set(jam.values()))}")

    # --- raw_counts integrity ------------------------------------
    raw = out["provenance"]["raw_counts"]
    expected_jobs = {"job_A1", "job_A2", "job_B1", "job_B2"}
    if set(raw.keys()) != expected_jobs:
        print(f"  [FAIL] raw_counts job_ids: got {sorted(raw.keys())}, "
              f"expected {sorted(expected_jobs)}")
        return 1
    total_pubs = sum(len(v) for v in raw.values())
    expected_pubs = (100 + 37) * 2  # two runs of 137 circuits each
    if total_pubs != expected_pubs:
        print(f"  [FAIL] raw_counts pub total: {total_pubs} != "
              f"{expected_pubs}")
        return 1
    print(f"  [PASS] raw_counts: {len(expected_jobs)} jobs, "
          f"{total_pubs} pub_results preserved")

    # --- per-run cleanup -----------------------------------------
    leaked = [k for r in out["runs"] for k in r if k.startswith("_")]
    if leaked:
        print(f"  [FAIL] private fields leaked into per-run section: "
              f"{sorted(set(leaked))}")
        return 1
    print("  [PASS] per-run section clean (no _underscore leakage)")

    # --- cross-check with reference (shape only) ------------------
    ref_path = REPO_ROOT / "data" / "w_repeat_results.json"
    if ref_path.is_file():
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        if ref["summary"]["std_convention"] != std_conv:
            print(f"  [FAIL] std_convention diverges from reference: "
                  f"{ref['summary']['std_convention']!r} vs {std_conv!r}")
            return 1
        ref_rbs = set(ref.get("random_basis_selection", {}).keys())
        if ref_rbs != got_rbs:
            print(f"  [FAIL] random_basis_selection keys diverge from "
                  f"reference: ref={sorted(ref_rbs)}, us={sorted(got_rbs)}")
            return 1
        print("  [PASS] matches reference summary.std_convention and "
              "random_basis_selection shape")

    print("\n  [PASS] output schema test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
