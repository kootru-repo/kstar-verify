"""Live-QPU replication of the flagship W-state K* measurement.

Submits 137 K*-selected Pauli circuits plus 137 seed-varying random
Pauli circuits per run against a user-chosen IBM Quantum backend,
then locally reconstructs both density matrices and reports
F(K*), F(random), and dF = F(K*) - F(random).

Requirements:
  - IBM Quantum Platform account with QPU time on a >= 4-qubit backend
  - IBM_QUANTUM_TOKEN environment variable set to the API token
  - Python deps:  pip install -r hardware/requirements.txt

Examples:
  # Dry run: print the protocol + QPU-time estimate, no submission.
  python hardware/run_w_state.py --dry-run

  # Minimal replication against the default least-busy backend.
  export IBM_QUANTUM_TOKEN=...
  python hardware/run_w_state.py --n-runs 2

  # Full flagship replication: 4 runs on ibm_fez (matches the shipped
  # reference at data/w_repeat_results.json).
  python hardware/run_w_state.py --backend ibm_fez --n-runs 4

  # Offline reanalysis of a saved JSON (no QPU time).
  python hardware/run_w_state.py --recover results.json

Exit codes:
  0  run completed (or --check-bounds PASS / --verify-equivalence match)
  1  --check-bounds FAIL, --verify-equivalence mismatch, or unrecoverable error
  2  bad CLI args, missing token env, missing cached files
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent

# Reuse the public repo's K* operator list and robust-MLE reconstruction
# rather than duplicating them inside hardware/.  The K* 137-operator
# set is shipped deterministically at data/kstar_operators_n4.json, and
# robust_mle / state_fidelity live in tier4-independent/.
sys.path.insert(0, str(REPO_ROOT / "tier4-independent"))


def _load_kstar_operators(n_qubits: int):
    """Return (pauli_ops_as_matrices, labels) for the n=4 K* set."""
    import core  # type: ignore  # from tier4-independent/core.py
    ops, labels, _ = core.select_kstar_paulis(n_qubits)
    return ops, labels


def _sample_random_operators(n_qubits: int, n_ops: int, seed: int):
    import core  # type: ignore
    ops, labels, _ = core.select_random_paulis(n_qubits, n_ops, seed=seed)
    return ops, labels


def _state_fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
    import core  # type: ignore
    return core.state_fidelity(rho, sigma)


def _reconstruct(expectations, pauli_ops, n_qubits, n_shots):
    import robust_mle  # type: ignore  # from tier4-independent/robust_mle.py
    return robust_mle.reconstruct_robust_mle(
        expectations, pauli_ops, n_qubits, n_shots=n_shots,
    )


def _default_seeds(n_runs: int):
    base = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    if n_runs > len(base):
        raise ValueError(f"--n-runs {n_runs} exceeds the built-in seed table of {len(base)}")
    return base[:n_runs]


def estimate_qpu_seconds(n_runs: int, n_shots: int, n_operators: int = 137) -> float:
    """Conservative estimate of billable QPU time in seconds.
    Assumes 2*n_operators circuits per run and ~1e-4 s per shot on a
    typical Eagle/Heron backend, plus a few seconds of per-run
    overhead.  Actual wall time is usually dominated by queue depth,
    not QPU time.
    """
    per_circuit = n_shots * 1e-4
    per_run = 2 * n_operators * per_circuit + 3.0
    return per_run * n_runs


def _connect(token: str, backend_name: str | None, n_qubits: int):
    from qiskit_ibm_runtime import QiskitRuntimeService

    service = QiskitRuntimeService(
        channel="ibm_quantum_platform",
        token=token,
    )
    if backend_name:
        return service, service.backend(backend_name)
    backends = service.backends(
        filters=lambda b: b.num_qubits >= n_qubits and b.status().operational,
    )
    if not backends:
        raise RuntimeError("No operational backend with enough qubits is available "
                           "on this account.  Pass --backend explicitly.")
    return service, min(backends, key=lambda b: b.status().pending_jobs)


def _account_identifier(service) -> str:
    """Return a non-sensitive account identifier for provenance.
    Falls back to a constant sentinel if qiskit_ibm_runtime does not
    expose an instance, so the reference JSON's key shape stays
    intact even on older runtimes / unusual auth setups.
    """
    try:
        acct = service.active_account() or {}
        # Never return the token.  Prefer instance (hub/group/project),
        # then channel, then a neutral sentinel.
        return (acct.get("instance")
                or acct.get("channel")
                or "active_account")
    except Exception:
        return "active_account"


def _capture_calibration(backend, physical_qubits: list[int]) -> dict | None:
    """Summarise backend.properties() for the physical qubits we
    actually use, matching the reference JSON's provenance.calibration
    shape (num_qubits / backend_name / n_qubit_properties /
    qubit_properties / n_gate_errors / gate_errors / retrieval_timestamp).

    Returns None if properties() is unavailable (FakeBackends,
    simulators, older runtimes) so _finalize can omit the field
    rather than emit a misleading empty object.
    """
    try:
        props = backend.properties()
        if props is None:
            return None
        d = props.to_dict()
    except Exception:
        return None

    qubits_list = d.get("qubits", []) or []
    gates_list = d.get("gates", []) or []

    # Slice per-qubit properties to just the physical qubits we used.
    used_qubits: list = []
    for q in physical_qubits:
        if 0 <= q < len(qubits_list):
            used_qubits.append(qubits_list[q])

    # Filter gates whose qubit support lies entirely inside our
    # physical_qubits set.  Two-qubit gates on neighbouring qubits
    # we use should be included; gates bridging to unused qubits
    # should not.
    phys_set = set(physical_qubits)
    used_gates = [
        g for g in gates_list
        if set(g.get("qubits", [])).issubset(phys_set)
    ]

    backend_name = (d.get("backend_name")
                    or getattr(backend, "name", "unknown"))
    return {
        "num_qubits": d.get("num_qubits") or len(qubits_list),
        "backend_name": backend_name,
        "n_qubit_properties": len(qubits_list),
        "qubit_properties": used_qubits,
        "n_gate_errors": len(gates_list),
        "gate_errors": used_gates,
        "retrieval_timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }


QPU_SECONDS_CEILING = 300.0  # 5 minutes


def _format_error(exc: Exception, context: str) -> tuple[str, str]:
    """Map a caught exception to (message, suggestion) strings.

    The classifier matches on exception-type name + message substring;
    it has no dependency on qiskit-ibm-runtime internals, so it keeps
    working across runtime version bumps (those tend to rename
    exception classes but keep the message wording recognisable).
    """
    name = type(exc).__name__
    msg = str(exc)
    low = msg.lower()

    # Authentication / authorization
    if ("auth" in low or "unauthorized" in low or "401" in low
            or "403" in low or "invalid token" in low
            or "not a valid api token" in low):
        return (
            f"Authentication failed during {context}: {name}: {msg}",
            "The API token was rejected by IBM Quantum.  Try:\n"
            "  1. Check the token has no leading/trailing whitespace:\n"
            "       echo \"[$IBM_QUANTUM_TOKEN]\"\n"
            "  2. Regenerate the token at https://quantum.ibm.com/account\n"
            "     (tokens expire; saved accounts go stale).\n"
            "  3. Confirm the token is bound to the ibm_quantum_platform\n"
            "     channel.  The older ibm_quantum channel was retired in 2024.\n"
            "  4. If you have multiple accounts, double-check which one\n"
            "     IBM_QUANTUM_TOKEN holds."
        )
    # Backend not found / no access
    if ("not found" in low or "no such backend" in low
            or "does not exist" in low or "backend.*unavailable" in low
            or "not accessible" in low):
        return (
            f"Backend not accessible during {context}: {name}: {msg}",
            "The requested backend is not reachable from your account.\n"
            "  1. List what IS available:\n"
            "       python -c \"from qiskit_ibm_runtime import "
            "QiskitRuntimeService; s = QiskitRuntimeService("
            "channel='ibm_quantum_platform'); "
            "print([b.name for b in s.backends()])\"\n"
            "  2. Drop --backend to let the runner auto-pick the "
            "least-busy 4+ qubit device.\n"
            "  3. Check your plan tier — some devices are paid-only."
        )
    # Network / timeout
    if (isinstance(exc, TimeoutError)
            or "timeout" in low or "timed out" in low
            or "connection" in low or "network" in low
            or "name or service not known" in low or "dns" in low):
        return (
            f"Network error during {context}: {name}: {msg}",
            "Could not reach IBM Quantum's API.\n"
            "  1. Check https://status.quantum.ibm.com for service outages.\n"
            "  2. Check your local network (VPN / proxy may block the API).\n"
            "  3. Retry — a transient timeout doesn't consume QPU budget."
        )
    # Quota / billing
    if ("quota" in low or "insufficient" in low or "not enough time" in low
            or "billing" in low):
        return (
            f"Quota or billing issue during {context}: {name}: {msg}",
            "Your free-tier QPU budget (10 min/month) may be exhausted.\n"
            "  1. Check usage at https://quantum.ibm.com/account.\n"
            "  2. Wait for the monthly reset, or upgrade your plan.\n"
            "  3. Cancel pending jobs that may be consuming the quota."
        )
    # Anything else
    return (
        f"Unexpected {name} during {context}: {msg}",
        "An unanticipated error was raised.  Diagnostic steps:\n"
        "  1. Rerun with --dry-run to confirm the protocol is well-formed.\n"
        "  2. Confirm dependency versions:\n"
        "       pip show qiskit qiskit-ibm-runtime qiskit-aer\n"
        "     (requirements.txt pins qiskit-ibm-runtime>=0.25,<2.0)\n"
        "  3. If the error persists, open an issue attaching the full\n"
        "     traceback and the 'pip show' output above."
    )


def _report_error(exc: Exception, context: str) -> int:
    """Print a formatted error + suggestion to stderr; return exit
    code 1.  All live-QPU code paths funnel through this helper so
    the user always sees both what broke and what to try next.
    """
    message, suggestion = _format_error(exc, context)
    print(f"\nERROR: {message}", file=sys.stderr)
    print(f"\nSUGGESTED FIX:\n{suggestion}", file=sys.stderr)
    return 1


def _tqdm(iterable=None, **kwargs):
    """Return a tqdm progress bar that works in terminal and notebook.
    Falls back to a no-op wrapper if tqdm is not installed (keeps the
    runner functional for an offline --dry-run / --recover path).
    """
    try:
        from tqdm.auto import tqdm
        return tqdm(iterable, **kwargs) if iterable is not None else tqdm(**kwargs)
    except ImportError:
        class _NoopBar:
            def __init__(self, it=None, **_): self._it = it
            def __iter__(self): return iter(self._it) if self._it is not None else iter([])
            def __enter__(self): return self
            def __exit__(self, *a): return False
            def update(self, *_): pass
            def set_postfix_str(self, *_): pass
            def close(self): pass
        return _NoopBar(iterable, **kwargs)


def _budget_gate(n_runs: int, n_shots: int, override: bool) -> int:
    """Refuse to start a run whose estimated QPU time exceeds
    5 minutes unless the caller explicitly acknowledged the budget
    with --yes-over-5min.  Protects the free-tier 10-min monthly
    allocation against accidental blow-outs from a mistyped flag.
    """
    est = estimate_qpu_seconds(n_runs, n_shots)
    if est <= QPU_SECONDS_CEILING or override:
        if override and est > QPU_SECONDS_CEILING:
            print(f"  [WARN] estimated QPU time {est:.0f}s > {QPU_SECONDS_CEILING:.0f}s; "
                  "proceeding because --yes-over-5min was passed.", flush=True)
        return 0
    print(
        f"ERROR: estimated QPU time is {est:.0f}s ({est/60:.1f} min), "
        f"which exceeds the {QPU_SECONDS_CEILING:.0f}s (5-min) safety gate.\n"
        "       This is protection against accidental budget blow-outs on the "
        "IBM free tier (10 min/month).\n"
        "       If you meant to do a longer run, re-run with --yes-over-5min.\n"
        "       Otherwise reduce --n-runs or --shots.",
        file=sys.stderr,
    )
    return 2


def _submit_and_read(backend, transpiled, n_shots: int, batch_size: int = 100):
    """Submit transpiled circuits in batches and return
    (batches, job_ids) where batches is a list of
    (batch_start_index, pub_result_list, job_id) tuples.  Pairing the
    job_id with its batch lets the caller build a provenance.raw_counts
    map keyed by job_id, matching the reference JSON layout.

    Raises the underlying exception on submit/result failure; callers
    wrap with _report_error so the user sees an actionable message
    instead of a bare traceback.
    """
    from qiskit_ibm_runtime import SamplerV2

    sampler = SamplerV2(mode=backend)
    batches = []
    job_ids = []

    n = len(transpiled)
    n_batches = (n + batch_size - 1) // batch_size
    with _tqdm(total=n_batches, desc="    submit batches",
               unit="batch", leave=False) as bar:
        for b, i in enumerate(range(0, n, batch_size), start=1):
            batch = transpiled[i:i + batch_size]
            pubs = [(circ, None, n_shots) for circ in batch]
            job = sampler.run(pubs)
            jid = job.job_id()
            job_ids.append(jid)
            bar.set_postfix_str(f"job={jid[:12]}")
            result = job.result()
            batches.append((i, result, jid))
            bar.update(1)
    return batches, job_ids


def _extract_counts(pub_result) -> dict:
    """Pull the measurement-register counts dict out of a PubResult,
    tolerant to the qiskit-runtime v2 data-bin layout.
    """
    try:
        return pub_result.data.meas.get_counts()
    except AttributeError:
        data_bin = pub_result.data
        for attr_name in dir(data_bin):
            attr = getattr(data_bin, attr_name)
            if hasattr(attr, "get_counts"):
                return attr.get_counts()
    raise RuntimeError("Could not locate a counts-producing field on pub_result.data")


def run_one(backend, kstar_circuits, kstar_ops, kstar_labels,
            rand_circuits, rand_ops, rand_labels,
            n_qubits: int, n_shots: int, rho_ideal: np.ndarray, run_idx: int, seed: int):
    """Submit one run (K* + random) and return a result dict."""
    from qiskit.transpiler import preset_passmanagers
    from circuits import counts_to_expectation

    all_circuits = kstar_circuits + rand_circuits
    print(f"  Circuits: {len(all_circuits)} "
          f"({len(kstar_circuits)} K* + {len(rand_circuits)} random)", flush=True)
    print("  Transpiling ...", flush=True)
    pm = preset_passmanagers.generate_preset_pass_manager(
        optimization_level=1, backend=backend,
    )
    transpiled = pm.run(all_circuits)
    # Capture the virtual-to-physical qubit mapping from the first
    # transpiled circuit (same layout applies to all 274).  Extract by
    # walking the measurement instructions rather than touching the
    # Layout API directly -- that keeps us robust across qiskit 1.x/2.x.
    physical_qubits: list[int] = []
    try:
        seen = set()
        for instr in transpiled[0].data:
            if instr.operation.name == "measure":
                phys = transpiled[0].find_bit(instr.qubits[0]).index
                if phys not in seen:
                    physical_qubits.append(phys)
                    seen.add(phys)
    except Exception:
        # Fallback: the runner will still emit -- just without this field.
        physical_qubits = []

    print("  Submitting ...", flush=True)
    batches, job_ids = _submit_and_read(backend, transpiled, n_shots)

    kstar_expectations: list[float] = []
    rand_expectations: list[float] = []
    # Per-job raw counts, matching the reference JSON's
    # provenance.raw_counts layout: {job_id: [counts_dict_per_pub]}.
    raw_counts_per_job: dict[str, list[dict]] = {}
    total_pubs = sum(len(r[1]) for r in batches)
    with _tqdm(total=total_pubs, desc="    extract counts",
               unit="circuit", leave=False) as bar:
        for batch_start, result, jid in batches:
            job_batch_counts: list[dict] = []
            for j in range(len(result)):
                idx = batch_start + j
                counts = _extract_counts(result[j])
                job_batch_counts.append(dict(counts))
                if idx < len(kstar_labels):
                    label = kstar_labels[idx]
                    kstar_expectations.append(
                        counts_to_expectation(counts, label, n_qubits))
                else:
                    label = rand_labels[idx - len(kstar_labels)]
                    rand_expectations.append(
                        counts_to_expectation(counts, label, n_qubits))
                bar.update(1)
            raw_counts_per_job[jid] = job_batch_counts

    print("  Reconstructing (local MLE) ...", flush=True)
    rho_kstar = _reconstruct(kstar_expectations, kstar_ops, n_qubits, n_shots)
    rho_rand = _reconstruct(rand_expectations, rand_ops, n_qubits, n_shots)

    f_kstar = float(_state_fidelity(rho_ideal, rho_kstar))
    f_rand = float(_state_fidelity(rho_ideal, rho_rand))
    delta_f = f_kstar - f_rand

    print(f"  F(K*)   = {f_kstar:.4f}", flush=True)
    print(f"  F(rand) = {f_rand:.4f}", flush=True)
    print(f"  dF      = {delta_f:+.4f}", flush=True)

    return {
        "run": run_idx + 1,
        "seed": seed,
        "f_kstar": f_kstar,
        "f_rand": f_rand,
        "delta_f": delta_f,
        "job_ids": job_ids,
        "kstar_expectations": [float(e) for e in kstar_expectations],
        "rand_expectations": [float(e) for e in rand_expectations],
        "rand_labels": rand_labels,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        # Private fields; lifted into the top-level `provenance` block
        # by _finalize so the output JSON layout matches the shipped
        # reference.  Underscore prefix marks them for strip-before-save.
        "_raw_counts_per_job": raw_counts_per_job,
        "_physical_qubits": physical_qubits,
    }


def run_w_repeat(token: str, backend_name: str | None, n_runs: int,
                 n_shots: int, output_dir: Path) -> dict:
    from circuits import pauli_label_to_circuit, ideal_state_for_prep

    n_qubits = 4
    seeds = _default_seeds(n_runs)
    kstar_ops, kstar_labels = _load_kstar_operators(n_qubits)
    rho_ideal = ideal_state_for_prep("w", n_qubits)

    print("=" * 70)
    print("  W-STATE K* — LIVE QPU REPLICATION")
    print(f"  runs = {n_runs}   shots = {n_shots}   operators = 2 x {len(kstar_labels)}")
    print(f"  seeds = {seeds}")
    print("=" * 70)
    print("  Connecting to IBM Quantum ...", flush=True)
    try:
        service, backend = _connect(token, backend_name, n_qubits)
    except Exception as exc:
        raise RuntimeError(
            _format_error(exc, "IBM Quantum connection")[0]
        ) from exc
    account_id = _account_identifier(service)
    print(f"  Backend: {backend.name}  (pending jobs: {backend.status().pending_jobs})",
          flush=True)

    print("  Building K* circuits (shared across runs) ...", flush=True)
    kstar_circuits = [
        pauli_label_to_circuit(label, n_qubits, prep="w") for label in kstar_labels
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    all_run_results: list[dict] = []
    calibration_snapshot: dict | None = None
    run_bar = _tqdm(total=n_runs, desc="runs", unit="run")
    try:
        for run_idx, seed in enumerate(seeds):
            print(f"\n{'=' * 70}\n  RUN {run_idx + 1}/{n_runs}   seed = {seed}\n{'=' * 70}",
                  flush=True)
            rand_ops, rand_labels = _sample_random_operators(n_qubits, len(kstar_ops), seed)
            rand_circuits = [
                pauli_label_to_circuit(label, n_qubits, prep="w") for label in rand_labels
            ]
            run_result = run_one(
                backend, kstar_circuits, kstar_ops, kstar_labels,
                rand_circuits, rand_ops, rand_labels,
                n_qubits=n_qubits, n_shots=n_shots, rho_ideal=rho_ideal,
                run_idx=run_idx, seed=seed,
            )
            # Capture backend calibration once, right after the first
            # successful transpile/submit pair -- physical_qubits is
            # populated on run_result by then, and subsequent runs share
            # the same calibration window (minutes apart).
            if calibration_snapshot is None:
                calibration_snapshot = _capture_calibration(
                    backend, run_result.get("_physical_qubits") or [],
                )
            # Propagate account_id via the private field so _finalize can
            # build provenance.job_account_map keyed by the per-run job_ids.
            run_result["_account_id"] = account_id
            all_run_results.append(run_result)
            _checkpoint(output_dir, backend.name, seeds, n_shots, n_runs, all_run_results)
            run_bar.update(1)
            run_bar.set_postfix_str(f"dF={run_result['delta_f']:+.3f}")
    finally:
        run_bar.close()

    return _finalize(output_dir, backend.name, seeds, n_shots, n_runs,
                     all_run_results, calibration=calibration_snapshot)


def _checkpoint(output_dir: Path, backend_name: str, seeds: list[int],
                n_shots: int, n_runs: int, runs: list[dict]) -> None:
    # Build a cumulative provenance.raw_counts map and a clean
    # per-run section *without* mutating `runs` itself -- the caller
    # keeps the original list for the next iteration, so we can't
    # pop() on the live dicts.  Copy shallowly and extract the
    # private field from each copy.
    raw_counts_all: dict[str, list[dict]] = {}
    runs_snapshot: list[dict] = []
    for r in runs:
        rr = dict(r)
        priv = rr.pop("_raw_counts_per_job", {})
        rr.pop("_physical_qubits", None)
        rr.pop("_account_id", None)
        raw_counts_all.update(priv)
        runs_snapshot.append(rr)

    payload = {
        "experiment": "w_state_repeat",
        "backend": backend_name,
        "n_qubits": 4,
        "n_shots": n_shots,
        "n_runs_completed": len(runs),
        "n_runs_planned": n_runs,
        "seeds": seeds,
        "runs": runs_snapshot,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "provenance": {"raw_counts": raw_counts_all},
    }
    # Include the completed-run count in the filename so two runs
    # finishing within the same clock-second don't overwrite each
    # other's checkpoint (unlikely on QPU but common in tests/sims).
    ck = output_dir / (
        f"w_repeat_partial_{time.strftime('%Y%m%d_%H%M%S')}_run{len(runs):02d}.json"
    )
    with ck.open("w", encoding="utf-8") as f:
        # default=str serialises datetime (and any other non-native
        # type backend.properties().to_dict() can embed inside
        # provenance.calibration) as an ISO string rather than
        # raising TypeError.
        json.dump(payload, f, indent=2, default=str)
    print(f"  [checkpoint] {ck.name}", flush=True)


def _finalize(output_dir: Path, backend_name: str, seeds: list[int],
              n_shots: int, n_runs: int, runs: list[dict],
              calibration: dict | None = None) -> dict:
    deltas = np.array([r["delta_f"] for r in runs])
    f_kstars = np.array([r["f_kstar"] for r in runs])
    f_rands = np.array([r["f_rand"] for r in runs])

    print(f"\n{'=' * 70}\n  SUMMARY ({n_runs} runs)\n{'=' * 70}")
    print(f"  {'run':>4s}  {'seed':>6s}  {'F(K*)':>8s}  {'F(rand)':>8s}  {'dF':>8s}")
    print(f"  {'-' * 4}  {'-' * 6}  {'-' * 8}  {'-' * 8}  {'-' * 8}")
    for r in runs:
        print(f"  {r['run']:>4d}  {r['seed']:>6d}  "
              f"{r['f_kstar']:>8.4f}  {r['f_rand']:>8.4f}  {r['delta_f']:>+8.4f}")
    # ddof=1 (sample std) matches the reference JSON's std_convention
    # field: "sample (ddof=1, N-1 denominator)".  On a single run, fall
    # back to ddof=0 to avoid division by zero -- an edge case the
    # reference never hits because it has 4 runs.
    ddof = 1 if len(runs) > 1 else 0
    print(f"\n  F(K*)   mean = {f_kstars.mean():.4f}  std = {f_kstars.std(ddof=ddof):.4f}")
    print(f"  F(rand) mean = {f_rands.mean():.4f}  std = {f_rands.std(ddof=ddof):.4f}")
    print(f"  dF      mean = {deltas.mean():+.4f}  std = {deltas.std(ddof=ddof):.4f}")

    summary = {
        "f_kstar_mean": float(f_kstars.mean()),
        "f_kstar_std":  float(f_kstars.std(ddof=ddof)),
        "f_rand_mean":  float(f_rands.mean()),
        "f_rand_std":   float(f_rands.std(ddof=ddof)),
        "delta_f_mean": float(deltas.mean()),
        "delta_f_std":  float(deltas.std(ddof=ddof)),
        # Verbatim string from the reference JSON so bit-exact diffs
        # against data/w_repeat_results.json's summary block match.
        "std_convention": ("sample (ddof=1, N-1 denominator)"
                           if ddof == 1 else "population (ddof=0, N denominator)"),
    }

    # Lift per-run private fields (_raw_counts_per_job, _physical_qubits)
    # into a single top-level provenance block keyed by job_id, matching
    # the reference JSON layout.  Strip underscore-prefixed keys from
    # each run dict on the way through so the per-run section stays
    # clean.
    raw_counts_all: dict[str, list[dict]] = {}
    all_job_ids: list[str] = []
    physical_qubits: list[int] = []
    job_account_map: dict[str, str] = {}
    runs_clean: list[dict] = []
    for r in runs:
        priv_counts = r.pop("_raw_counts_per_job", {})
        priv_pq = r.pop("_physical_qubits", [])
        priv_acct = r.pop("_account_id", "active_account")
        raw_counts_all.update(priv_counts)
        all_job_ids.extend(r.get("job_ids", []))
        if priv_pq and not physical_qubits:
            # Same layout is used across runs (same backend, same
            # transpile args), so the first non-empty value wins.
            physical_qubits = list(priv_pq)
        for jid in r.get("job_ids", []):
            job_account_map[jid] = priv_acct
        runs_clean.append(r)

    provenance = {
        "raw_counts": raw_counts_all,
        "job_ids": all_job_ids,
        "n_jobs_retrieved": len(all_job_ids),
        "n_jobs_failed": 0,
        "job_account_map": job_account_map,
        "transpilation": {
            "optimization_level": 1,
            "note": ("All circuits transpiled with optimization_level=1 via "
                     "qiskit.transpiler.generate_preset_pass_manager."),
        },
        "retrieval_timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if physical_qubits:
        provenance["physical_qubits"] = physical_qubits
        provenance["physical_qubits_note"] = (
            "Virtual-to-physical qubit mapping chosen by the transpile "
            "preset pass manager at optimization_level=1; extracted "
            "from the first measurement instruction in each circuit."
        )
    if calibration is not None:
        provenance["calibration"] = calibration

    output = {
        "experiment": "w_state_repeat",
        "backend": backend_name,
        "n_qubits": 4,
        "n_shots": n_shots,
        "n_runs": n_runs,
        "seeds": seeds,
        "runs": runs_clean,
        "summary": summary,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        # Top-level metadata fields matching data/w_repeat_results.json
        # so a side-by-side JSON diff of a rerun on the same backend
        # yields only calibration-driven numeric differences, not
        # schema drift.
        "provenance": provenance,
        "random_basis_selection": {
            "method": ("select_random_paulis(n_qubits=4, "
                       "n_measurements=137, seed=<run_seed>)"),
            "per_run_seeds": seeds,
            "note": (
                "Each run uses its own seed for independent random "
                "operator selection.  The same 137 K* operators are "
                "used across all runs (deterministic)."
            ),
            "code": ("from core import select_random_paulis; "
                     "ops, labels, _ = select_random_paulis(4, 137, seed=seed)"),
        },
    }
    outfile = output_dir / f"w_repeat_results_{time.strftime('%Y%m%d_%H%M%S')}.json"
    with outfile.open("w", encoding="utf-8") as f:
        # default=str handles the datetime values that
        # backend.properties().to_dict() embeds inside
        # provenance.calibration.gate_errors / qubit_properties.
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Saved: {outfile}")
    return output


def reanalyze(path: Path) -> int:
    """Recompute summary statistics + bootstrap CI from a saved JSON."""
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    runs = data.get("runs", [])
    if not runs:
        print(f"  No runs found in {path}")
        return 2
    deltas = np.array([r["delta_f"] for r in runs])
    print(f"\n  Loaded {len(runs)} runs from {path}")
    print(f"  {'run':>4s}  {'seed':>6s}  {'F(K*)':>8s}  {'F(rand)':>8s}  {'dF':>8s}")
    print(f"  {'-' * 4}  {'-' * 6}  {'-' * 8}  {'-' * 8}  {'-' * 8}")
    for r in runs:
        print(f"  {r['run']:>4d}  {r['seed']:>6d}  "
              f"{r['f_kstar']:>8.4f}  {r['f_rand']:>8.4f}  {r['delta_f']:>+8.4f}")
    ddof = 1 if len(runs) > 1 else 0
    print(f"\n  dF mean = {deltas.mean():+.4f}   std = {deltas.std(ddof=ddof):.4f}")

    # Bootstrap 95% CI on dF mean.  Seeded for reproducibility.
    rng = np.random.default_rng(0)
    boot = rng.choice(deltas, size=(10_000, len(deltas)), replace=True).mean(axis=1)
    print(f"  95% CI = [{np.percentile(boot, 2.5):+.4f}, "
          f"{np.percentile(boot, 97.5):+.4f}]")
    print(f"  runs with dF > 0: {int((deltas > 0).sum())}/{len(deltas)}")
    return 0


def dry_run(n_runs: int, n_shots: int) -> int:
    est = estimate_qpu_seconds(n_runs, n_shots)
    print("DRY RUN — no QPU submission.")
    print(f"  Planned: {n_runs} run(s) x 274 circuits x {n_shots} shots")
    print(f"  Rough QPU-time estimate: {est:.1f} s ({est / 60:.2f} min) "
          "[excludes queue wait]")
    print("  Per run, the script does:")
    print("    1. Generate 137 K* Pauli labels deterministically via")
    print("       tier4-independent.core.select_kstar_paulis(4)")
    print("    2. Sample 137 random Pauli labels (seed-dependent)")
    print("    3. Build+transpile 274 circuits for the chosen backend")
    print("    4. Submit in batches of 100; record raw counts + job IDs")
    print("    5. Locally reconstruct rho_K* and rho_rand via robust MLE")
    print("    6. Report F(K*), F(rand), dF = F(K*) - F(rand)")
    print("    7. Checkpoint to JSON after every run")
    # Load the acceptance bounds and echo them so the dry-run reflects
    # whatever is actually in expected_bounds.json (single source of
    # truth) rather than hard-coded values that can drift apart.
    bounds_path = THIS_DIR / "expected_bounds.json"
    if bounds_path.is_file():
        with bounds_path.open(encoding="utf-8") as f:
            b = json.load(f)
        print()
        print("  Acceptance bounds (hardware/expected_bounds.json):")
        print(f"    F(K*) mean >= {b.get('f_kstar_min', 'n/a')}")
        print(f"    dF mean    >= {b.get('delta_f_min', 'n/a')}")
        print(f"    runs with dF > 0 >= {b.get('n_positive_runs_min', 'n/a')}")
    return 0


def parse_args(argv):
    p = argparse.ArgumentParser(
        description=("Live-QPU replication of the flagship 4-qubit W-state K* "
                     "measurement.  Reads IBM_QUANTUM_TOKEN from the environment."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--backend", default=None,
                   help="IBM backend name (default: least-busy >= 4 qubits)")
    p.add_argument("--n-runs", type=int, default=5,
                   help="Number of independent runs (default 5, matching the "
                        "original ibm_w_repeat.py default; pass --n-runs 4 to "
                        "match the run-count in data/w_repeat_results.json "
                        "exactly).")
    p.add_argument("--shots", type=int, default=1000,
                   help="Shots per circuit (default 1000)")
    p.add_argument("--output-dir", default=str(THIS_DIR / "results"),
                   help="Directory to write result JSON into")
    p.add_argument("--dry-run", action="store_true",
                   help="Print protocol and QPU estimate only, do not submit")
    p.add_argument("--recover", metavar="FILE",
                   help="Reanalyze a saved results JSON (no QPU submission)")
    p.add_argument("--verify-equivalence", action="store_true",
                   help="Check that core.select_kstar_paulis(4) matches the "
                        "cached data/kstar_operators_n4.json byte-for-byte; "
                        "no QPU submission")
    p.add_argument("--preflight", action="store_true",
                   help="Submit a 2-circuit sanity check (~2 s QPU time) "
                        "that validates token + backend + SamplerV2 + "
                        "counts extraction + basis rotation before the "
                        "user commits the full ~2.5-min run.  Exits after "
                        "the check; re-invoke without --preflight for the "
                        "full replication.")
    p.add_argument("--check-bounds", action="store_true",
                   help="After a completed run, evaluate expected_bounds.json "
                        "and exit non-zero on FAIL.  Default behavior is to "
                        "save the JSON and exit 0 regardless of numbers -- "
                        "the bounds layer is an interpretation not part of "
                        "the raw protocol.")
    p.add_argument("--yes-over-5min", action="store_true",
                   help="Acknowledge that the estimated QPU time exceeds the "
                        "5-minute safety gate.  Required whenever --n-runs x "
                        "--shots would push the estimate above 300 s; "
                        "protects the IBM free-tier 10-min monthly allocation "
                        "against accidental blow-outs from a mistyped flag.")
    p.add_argument("--token-env", default="IBM_QUANTUM_TOKEN",
                   help="Env var holding the API token (default IBM_QUANTUM_TOKEN)")
    return p.parse_args(argv)


PREFLIGHT_CHECKS = [
    # (pauli label, noiseless <P>_|W_4>, tolerance)
    # <IIIZ>  = <Z_0>  on |W_4> = 0.5   (q0 is |1> in 1 of 4 basis states)
    # <IIXI>  = <X_1>  on |W_4> = 0.0   (W is orthogonal to its own X-flip)
    # Hardware noise lowers |<Z>| and keeps |<X>| ~0; tolerance 0.3 catches
    # gross pipeline breakage without false-failing on a noisy backend.
    ("IIIZ", 0.5, 0.3),
    ("IIXI", 0.0, 0.3),
]


def preflight(token: str, backend_name: str | None = None,
              n_shots: int = 1000) -> int:
    """Submit a 2-circuit sanity check that validates the full
    end-to-end pipeline before the user commits a full run's QPU
    budget.  At n_shots=1000 this costs roughly 0.2 s of actual QPU
    time plus per-submission overhead, well under the 5-minute gate.

    Tests that matter for a first-time user:
      - token authenticates against qiskit_ibm_runtime
      - chosen backend is reachable
      - SamplerV2(mode=backend) round-trips a pub successfully
      - pub_result.data.meas.get_counts() returns a dict
      - counts_to_expectation lands within a physically-reasonable
        band for the W-state prep (0.2-0.8 for <Z_0>, |v| < 0.3 for
        <X_1>)

    Returns 0 on all checks PASS, 1 on any FAIL.  Prints per-check
    detail so the user can see what specifically failed.
    """
    n_qubits = 4
    print("=" * 70)
    print("  QPU PREFLIGHT — 2-circuit sanity check (~2 s QPU budget)")
    print("=" * 70)

    from circuits import pauli_label_to_circuit, counts_to_expectation
    from qiskit.transpiler import preset_passmanagers
    from qiskit_ibm_runtime import SamplerV2

    # Phase 1: connect -------------------------------------------------
    print("  [1/4] Connecting to IBM Quantum ...", flush=True)
    try:
        _, backend = _connect(token, backend_name, n_qubits)
    except Exception as exc:
        return _report_error(exc, "IBM Quantum connection")
    print(f"        backend: {backend.name}  "
          f"(pending jobs: {backend.status().pending_jobs})", flush=True)

    # Phase 2: build + transpile ---------------------------------------
    print(f"  [2/4] Building and transpiling {len(PREFLIGHT_CHECKS)} "
          "circuits ...", flush=True)
    try:
        labels = [label for label, _, _ in PREFLIGHT_CHECKS]
        circuits = [pauli_label_to_circuit(lbl, n_qubits, prep="w") for lbl in labels]
        pm = preset_passmanagers.generate_preset_pass_manager(
            optimization_level=1, backend=backend,
        )
        transpiled = pm.run(circuits)
    except Exception as exc:
        return _report_error(exc, "circuit build / transpile")

    # Phase 3: submit + wait for result --------------------------------
    print(f"  [3/4] Submitting {len(circuits)} x {n_shots} shots "
          "(~2 s of QPU, plus queue wait) ...", flush=True)
    try:
        sampler = SamplerV2(mode=backend)
        pubs = [(c, None, n_shots) for c in transpiled]
        with _tqdm(total=1, desc="        submit", unit="job", leave=False) as bar:
            job = sampler.run(pubs)
            jid = job.job_id()
            bar.set_postfix_str(f"job={jid[:12]}")
            result = job.result()
            bar.update(1)
        print(f"        job id: {jid}", flush=True)
    except Exception as exc:
        return _report_error(exc, "sampler submission")

    # Phase 4: extract + evaluate --------------------------------------
    print(f"  [4/4] Extracting counts and computing expectations ...",
          flush=True)
    try:
        all_pass = True
        rows: list[tuple[str, bool, float, float, float]] = []
        with _tqdm(total=len(PREFLIGHT_CHECKS), desc="        extract",
                   unit="circuit", leave=False) as bar:
            for idx, (label, expected, tol) in enumerate(PREFLIGHT_CHECKS):
                counts = _extract_counts(result[idx])
                measured = counts_to_expectation(counts, label, n_qubits)
                ok = abs(measured - expected) < tol
                all_pass &= ok
                rows.append((label, ok, measured, expected, tol))
                bar.update(1)
    except Exception as exc:
        return _report_error(exc, "counts extraction")

    # Report ---------------------------------------------------------
    print()
    for label, ok, measured, expected, tol in rows:
        tag = "PASS" if ok else "FAIL"
        print(f"  [{tag}] <{label}>  measured = {measured:+.4f}  "
              f"(noiseless = {expected:+.2f}, tol = {tol})")

    print()
    if all_pass:
        print("  [PASS] preflight succeeded -- pipeline wired correctly")
        print("         next: drop --preflight to submit the full run")
        return 0
    print("  [FAIL] preflight failed -- measured values outside tolerance")
    print()
    print("         Likely causes, in rough order of probability:")
    print("           a. Backend calibration has drifted since the "
          "reference run.  Re-check the IBM dashboard for recent")
    print("              calibration timestamps on the selected backend.")
    print("           b. Transpile picked different physical qubits than "
          "the [0,1,2,3] used for the shipped reference.  The run")
    print("              still works, but absolute F(K*) / dF numbers "
          "will shift.  To force specific qubits, pre-transpile with")
    print("              initial_layout=[0,1,2,3] or pick a backend where "
          "the first 4 qubits are the healthiest.")
    print("           c. tier4-independent/core.py's counts_to_expectation "
          "was changed and no longer matches this runner.  Run")
    print("              `python hardware/run_w_state.py --verify-equivalence`")
    print("              to confirm the K* operator set is still in sync.")
    return 1


def verify_equivalence() -> int:
    """Compare core.select_kstar_paulis(4) output against the cached
    export at data/kstar_operators_n4.json.  Both are expected to
    produce the same 137 Pauli labels in the same order; a mismatch
    means either the K* algorithm has drifted from what was shipped
    (code change without regenerating the JSON) or the JSON was
    replaced with a different set.  Returns 0 on match, 1 on drift.
    """
    import core  # type: ignore

    cached_path = REPO_ROOT / "data" / "kstar_operators_n4.json"
    if not cached_path.is_file():
        print(f"ERROR: cached K* export not found at {cached_path}",
              file=sys.stderr)
        return 2
    with cached_path.open(encoding="utf-8") as f:
        cached = json.load(f)
    cached_labels: list[str] = cached.get("operators", [])

    _, live_labels, _ = core.select_kstar_paulis(4)

    print("Protocol equivalence check")
    print(f"  cached file: {cached_path}")
    print(f"  cached count: {len(cached_labels)}")
    print(f"  live count:   {len(live_labels)}")

    if cached_labels == live_labels:
        print("  [PASS] 137 K* Pauli labels match in identical order.")
        return 0

    print("  [FAIL] K* label sequence differs between code and cache.")
    # Surface the first few differences for quick debugging.
    shown = 0
    for i, (a, b) in enumerate(zip(cached_labels, live_labels)):
        if a != b:
            print(f"    index {i:3d}: cached={a!r}  live={b!r}")
            shown += 1
            if shown >= 10:
                break
    if len(cached_labels) != len(live_labels):
        print(f"    length differs by {abs(len(cached_labels) - len(live_labels))}")
    return 1


def main(argv=None):
    args = parse_args(argv)
    if args.verify_equivalence:
        return verify_equivalence()
    if args.recover:
        return reanalyze(Path(args.recover))
    if args.dry_run:
        return dry_run(args.n_runs, args.shots)

    token = os.environ.get(args.token_env)
    if not token:
        print(f"ERROR: env var {args.token_env} is unset.  "
              f"Get a token at https://quantum.ibm.com and run:", file=sys.stderr)
        print(f"  export {args.token_env}=<your-token>", file=sys.stderr)
        return 2

    if args.preflight:
        return preflight(token, args.backend)

    # Budget gate BEFORE we connect so a user who'd blow their free tier
    # hears about it immediately, not after queueing jobs.
    gate_rc = _budget_gate(args.n_runs, args.shots, args.yes_over_5min)
    if gate_rc != 0:
        return gate_rc

    try:
        output = run_w_repeat(
            token=token,
            backend_name=args.backend,
            n_runs=args.n_runs,
            n_shots=args.shots,
            output_dir=Path(args.output_dir),
        )
    except Exception as exc:
        # run_w_repeat rewrapped its own errors as RuntimeError with the
        # _format_error message already.  Everything else is an
        # unexpected exception -- surface both with the same handler.
        return _report_error(exc, "live-run execution")
    if args.check_bounds:
        return _check_all_bounds(output)
    # Default path: emit the raw numbers and let the user (or
    # analyze_results.py) apply whatever acceptance criteria they
    # want.  Exit 0 regardless of the numbers so the run is treated
    # as a successful submission + JSON capture.
    return 0


def _check_all_bounds(output: dict) -> int:
    """Evaluate every threshold in expected_bounds.json and return
    0 on all-PASS, 1 on any-FAIL.  Kept separate from analyze_results.py
    so the runner stays importable without the analyzer, but the two
    must stay in sync (same bounds JSON, same three checks).
    """
    bounds_path = THIS_DIR / "expected_bounds.json"
    if not bounds_path.is_file():
        return 0
    with bounds_path.open(encoding="utf-8") as f:
        b = json.load(f)
    summary = output["summary"]
    n_positive = sum(1 for r in output["runs"] if r["delta_f"] > 0)
    checks = [
        ("F(K*) mean",       summary["f_kstar_mean"] >= b.get("f_kstar_min", 0.0),
         summary["f_kstar_mean"], b.get("f_kstar_min", 0.0)),
        ("dF mean",          summary["delta_f_mean"] >= b.get("delta_f_min", 0.0),
         summary["delta_f_mean"], b.get("delta_f_min", 0.0)),
        ("runs with dF > 0", n_positive >= b.get("n_positive_runs_min", 0),
         n_positive, b.get("n_positive_runs_min", 0)),
    ]
    print("\n  Bound check (hardware/expected_bounds.json):")
    all_ok = True
    for name, ok, got, need in checks:
        all_ok &= ok
        print(f"    [{'PASS' if ok else 'FAIL'}] {name}: got {got}, need >= {need}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
