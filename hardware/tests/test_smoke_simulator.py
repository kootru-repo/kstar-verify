"""Local-simulator smoke test for hardware/run_w_state.py.

Exercises every qiskit code path the live-QPU runner uses, but
against a local simulator backend (FakeBrisbane via qiskit-ibm-
runtime's fake_provider) so no QPU credits or tokens are spent:

  hardware/circuits.py  -> pauli_label_to_circuit
                        -> counts_to_expectation
                        -> ideal_state_for_prep
  qiskit.transpiler     -> generate_preset_pass_manager
  qiskit_ibm_runtime    -> SamplerV2(mode=fake_backend)
                        -> job.result() -> pub_result.data.meas
  tier4-independent     -> robust_mle.reconstruct_robust_mle
                        -> core.state_fidelity
                        -> core.select_kstar_paulis
                        -> core.select_random_paulis

A smaller subset (20 + 20 operators, 500 shots) is used for speed.
FakeBrisbane adds realistic noise, so F(K*) will NOT be 1.0; the
assertion is only that the pipeline produces plausible, finite
numbers and that F(K*) > F(random) on this noiseless-plus-device-
noise simulator.  A stricter end-to-end check would swap FakeBrisbane
for qiskit_aer.AerSimulator, but the fake provider is what validates
the exact SamplerV2(mode=backend) code path the real runner uses.

Skipped if qiskit / qiskit-ibm-runtime are not installed.

Usage:
  python hardware/tests/test_smoke_simulator.py
"""
from __future__ import annotations

import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
HARDWARE_DIR = THIS_DIR.parent
REPO_ROOT = HARDWARE_DIR.parent

# Make circuits.py, tier4 core.py, tier4 robust_mle.py importable.
for p in (HARDWARE_DIR, REPO_ROOT / "tier4-independent"):
    sys.path.insert(0, str(p))


def _skip(msg: str) -> int:
    print(f"SKIP: {msg}")
    return 0


def main() -> int:
    try:
        from qiskit_ibm_runtime import SamplerV2
        from qiskit_ibm_runtime.fake_provider import FakeBrisbane
        from qiskit.transpiler import generate_preset_pass_manager
    except ImportError as e:
        return _skip(f"qiskit stack not installed ({e}); "
                     "run `pip install -r hardware/requirements.txt` first.")

    import numpy as np

    # The helpers under test.
    from circuits import (  # type: ignore  # under HARDWARE_DIR
        pauli_label_to_circuit,
        counts_to_expectation,
        ideal_state_for_prep,
    )
    import core  # type: ignore  # under tier4-independent
    import robust_mle  # type: ignore  # under tier4-independent

    n_qubits = 4
    n_shots = 500
    # Full 137-operator set: 20-of-137 leaves the MLE badly
    # underdetermined and F(K*) floors near 0.2 regardless of pipeline
    # correctness, making assertions noisy.  Full set takes ~30 s on a
    # local FakeBrisbane and yields ~0.7-0.85 F(K*) -- a meaningful
    # signal.
    n_subset: int | None = None  # None = full set

    print("=" * 70)
    print("  SIMULATOR SMOKE TEST — hardware/ end-to-end (no QPU credits)")
    print("  backend  = FakeBrisbane (local, noisy)")
    print(f"  operators = {'full 137' if n_subset is None else n_subset} K* + same count random")
    print(f"  shots    = {n_shots}")
    print("=" * 70)

    # --- Build inputs exactly as the runner does. ---
    kstar_ops_all, kstar_labels_all, _ = core.select_kstar_paulis(n_qubits)
    if n_subset is None:
        kstar_ops = kstar_ops_all
        kstar_labels = kstar_labels_all
    else:
        kstar_ops = kstar_ops_all[:n_subset]
        kstar_labels = kstar_labels_all[:n_subset]
    rand_ops, rand_labels, _ = core.select_random_paulis(
        n_qubits, len(kstar_ops), seed=100,
    )
    rho_ideal = ideal_state_for_prep("w", n_qubits)

    kstar_circuits = [
        pauli_label_to_circuit(lbl, n_qubits, prep="w") for lbl in kstar_labels
    ]
    rand_circuits = [
        pauli_label_to_circuit(lbl, n_qubits, prep="w") for lbl in rand_labels
    ]
    all_circuits = kstar_circuits + rand_circuits
    print(f"  [OK] built {len(all_circuits)} QuantumCircuit objects")

    # --- Transpile against a real backend model. ---
    backend = FakeBrisbane()
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
    transpiled = pm.run(all_circuits)
    print(f"  [OK] transpiled for {backend.name} ({backend.num_qubits} qubits)")

    # --- Submit via SamplerV2(mode=fake_backend): mirrors the real runner. ---
    sampler = SamplerV2(mode=backend)
    pubs = [(c, None, n_shots) for c in transpiled]
    job = sampler.run(pubs)
    result = job.result()
    print(f"  [OK] sampler.run returned {len(result)} pub results")

    # --- Extract counts and expectations (same logic as runner). ---
    kstar_exps: list[float] = []
    rand_exps: list[float] = []
    split = len(kstar_labels)
    for idx in range(len(result)):
        pub = result[idx]
        try:
            counts = pub.data.meas.get_counts()
        except AttributeError:
            data_bin = pub.data
            counts = None
            for attr_name in dir(data_bin):
                attr = getattr(data_bin, attr_name)
                if hasattr(attr, "get_counts"):
                    counts = attr.get_counts()
                    break
            if counts is None:
                raise RuntimeError(
                    f"could not extract counts from pub[{idx}]")
        if idx < split:
            label = kstar_labels[idx]
            kstar_exps.append(counts_to_expectation(counts, label, n_qubits))
        else:
            label = rand_labels[idx - split]
            rand_exps.append(counts_to_expectation(counts, label, n_qubits))
    print(f"  [OK] extracted {len(kstar_exps)} K* + {len(rand_exps)} random expectations")

    # --- Local reconstruction + fidelity (tier4). ---
    rho_kstar = robust_mle.reconstruct_robust_mle(
        kstar_exps, kstar_ops, n_qubits, n_shots=n_shots,
    )
    rho_rand = robust_mle.reconstruct_robust_mle(
        rand_exps, rand_ops, n_qubits, n_shots=n_shots,
    )
    f_kstar = float(core.state_fidelity(rho_ideal, rho_kstar))
    f_rand = float(core.state_fidelity(rho_ideal, rho_rand))
    delta_f = f_kstar - f_rand
    print("\n  F(K*)   = {:.4f}".format(f_kstar))
    print("  F(rand) = {:.4f}".format(f_rand))
    print("  dF      = {:+.4f}".format(delta_f))

    # --- Sanity assertions. ---
    # The smoke test validates the pipeline, not the flagship numbers.
    # FakeBrisbane's device noise typically drops F(K*) to ~0.7-0.85
    # on the full 137-op set, but F(K*) must still beat F(random) --
    # that's the one invariant that only holds if every link in the
    # pipeline (circuits, transpile, sampler, counts, MLE, fidelity)
    # is correct.
    assert np.isfinite(f_kstar) and 0.0 <= f_kstar <= 1.0, \
        f"F(K*) out of range: {f_kstar}"
    assert np.isfinite(f_rand) and 0.0 <= f_rand <= 1.0, \
        f"F(rand) out of range: {f_rand}"
    assert delta_f > 0.0, (
        f"dF={delta_f:+.3f} is non-positive: K* should beat random on a "
        "noiseless-plus-device-noise simulator.  Likely pipeline bug."
    )
    print("\n  [PASS] pipeline end-to-end validated on local simulator")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
