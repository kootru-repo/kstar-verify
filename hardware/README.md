# Live-QPU replication — flagship 4-qubit W-state K* measurement

This directory lets anyone with an IBM Quantum account and a few
minutes of QPU time reproduce the flagship hardware finding from the
companion manuscript:

> On a 4-qubit W state, the 137-element K*-selected measurement set
> recovers F(K\*) = 0.872 ± 0.021, whereas a same-size random Pauli set
> recovers F(rand) = 0.542 ± 0.076, giving an advantage
> ΔF = +0.330 ± 0.071 across 4 independent runs on `ibm_fez`.

The reference dataset (the JSON this protocol produced in March 2026)
ships at [`../data/w_repeat_results.json`](../data/w_repeat_results.json).
Running this module on your own account should land in the same
neighbourhood modulo backend calibration differences.

## What this runs

Per run (default 4 runs):
1. Generates the deterministic 137-operator K\* set via
   `core.select_kstar_paulis(4)` from
   [`../tier4-independent/core.py`](../tier4-independent/core.py).
   (The set is also cached at
   [`../data/kstar_operators_n4.json`](../data/kstar_operators_n4.json)
   for offline inspection; the runner recomputes it so a corrupted
   data file cannot silently change the protocol.)
2. Samples 137 random Pauli operators at a fixed seed
   (seeds 100, 200, 300, 400 by default; override via `--n-runs`).
3. Builds `2 × 137 = 274` measurement circuits on top of a
   hard-wired 4-qubit W-state preparation.
4. Transpiles + submits to the chosen IBM backend in batches of 100.
5. Locally reconstructs ρ\_K\* and ρ\_rand via robust MLE
   (reusing [`../tier4-independent/robust_mle.py`](../tier4-independent/robust_mle.py)).
6. Reports `F(K*)`, `F(rand)`, `ΔF = F(K*) − F(rand)`.
7. Checkpoints a JSON file after every completed run so a dropped
   connection does not waste earlier QPU time.

The protocol never sends anything beyond the 274 measurement circuits
and their shots — no external model, no hybrid optimization, no
mid-circuit feedback.

## Protocol equivalence with the flagship run

The runner is designed so that a fresh submission against the same
backend at the same calibration epoch would produce a JSON
statistically indistinguishable from
[`../data/w_repeat_results.json`](../data/w_repeat_results.json).
Everything that could change a number is held fixed:

| element | how it is held fixed |
| --- | --- |
| K\* operator set | `core.select_kstar_paulis(4)` — deterministic; cached at `../data/kstar_operators_n4.json` |
| random Pauli seeds | `[100, 200, 300, 400, 500]` (default `--n-runs 5`, matching the original script); pass `--n-runs 4` to match the four runs recorded in the shipped reference exactly |
| random operator sampling | `core.select_random_paulis(4, 137, seed=<run_seed>)` |
| W-state preparation | hard-wired CRy/CNOT cascade with splitting ratios 3/4, 2/3, 1/2 |
| measurement basis | per-qubit H / S†·H / identity for X / Y / Z |
| shots | `--shots 1000` default |
| transpile level | `optimization_level=1` via `generate_preset_pass_manager` |
| submission | `SamplerV2(mode=backend)`, batch size 100 |
| reconstruction | `tier4-independent/robust_mle.reconstruct_robust_mle` (same code as tier 4 of the verification suite) |
| fidelity | `core.state_fidelity` (Uhlmann fidelity) |
| summary std | `ddof=1` (sample std, `"sample (ddof=1, N-1 denominator)"`), matches the reference's `std_convention` |
| output schema | emits all 10 reference sub-keys under `provenance` (`raw_counts`, `job_ids`, `n_jobs_retrieved`, `n_jobs_failed`, `job_account_map`, `transpilation`, `retrieval_timestamp`, `physical_qubits`, `physical_qubits_note`, `calibration`) plus the 4-key `random_basis_selection` block — full structural parity |

Confirm end-to-end equivalence in one command:

```bash
python hardware/run_w_state.py --verify-equivalence
# expected: "[PASS] 137 K* Pauli labels match in identical order."
```

The only deliberate departures from the original private script are
packaging choices that cannot affect results:

1. **Token** — read from `IBM_QUANTUM_TOKEN` env var instead of
   `sys.argv[1]`.  Prevents the token from landing in shell history
   or process listings.
2. **Classical register** — named `"meas"` instead of the default
   `"c"` so `pub_result.data.meas.get_counts()` works on the happy
   path.  Register names are metadata, not quantum state; counts are
   bit-identical.
3. **Output directory** — defaults to `hardware/results/` instead of
   the current working directory; configurable via `--output-dir`.
4. **Checkpoint filename** — includes the run index
   (`_run01`, `_run02`, ...) so two runs finishing in the same
   clock-second cannot overwrite each other.
5. **Bound check** — opt-in via `--check-bounds`.  Default exit is 0
   regardless of the numbers; the user can apply their own
   acceptance criteria without the protocol's exit code editorializing.

For a byte-for-byte diff, run `analyze_results.py --reference` on the
shipped JSON and `analyze_results.py results/w_repeat_results_*.json`
on your own run — the fields printed are identical.

## Budget safety

The runner enforces a **5-minute QPU-time ceiling** per invocation
(`QPU_SECONDS_CEILING = 300.0` in `run_w_state.py`).  If
`--n-runs × --shots` would push the estimated QPU time above 300 s,
the script aborts with a clear error before connecting, so a mistyped
flag (`--n-runs 40`) cannot silently burn the free tier's 10-min
monthly allocation.  Override with `--yes-over-5min` when you
genuinely want a longer run.

At defaults (`--n-runs 5 --shots 1000`) the estimate is ~2.5 min —
safely under the gate and roughly one-quarter of the free-tier budget,
so the same account can do ~4 clean runs per month.

## Prerequisites

1. **IBM Quantum Platform account** — sign up at
   <https://quantum.ibm.com>.  The free open plan includes 10 minutes
   per month of QPU time on the Open backends, which is enough for
   ≥ 5 runs at 1000 shots (default configuration uses ~2.5 min).
2. **API token** — copy the token from the account dashboard and
   export it.  The runner reads the `IBM_QUANTUM_TOKEN` env var by
   default:
   ```bash
   export IBM_QUANTUM_TOKEN='your-token-here'
   ```
   The token itself is never written to disk; pass a different env
   var name via `--token-env NAME` if you prefer.
3. **Python deps** (into the same venv used by the rest of this repo):
   ```bash
   pip install -r requirements.txt              # top-level repo deps
   pip install -r hardware/requirements.txt     # adds qiskit-ibm-runtime
   ```

## One-command replication

From the repository root:

```bash
# Dry run first: prints protocol summary + QPU-time estimate,
# submits nothing.
python hardware/run_w_state.py --dry-run

# Preflight: 2 circuits, ~2 s of real QPU, validates token +
# backend + full pipeline before you commit a full run's budget.
python hardware/run_w_state.py --preflight

# Full flagship replication (matches the shipped reference's run count):
python hardware/run_w_state.py --backend ibm_fez --n-runs 4

# Shorter run on the least-busy available backend:
python hardware/run_w_state.py --n-runs 2
```

**Recommended first-run order:** `--dry-run` → `--preflight` → full run.
The preflight fails in seconds (not minutes) if anything is wrong,
saving the rest of your free-tier budget for the actual replication.

The same preflight is also wired into the Jupyter notebook as an
optional appendix (last two cells of `notebook/k_star_demo.ipynb`);
that appendix skips cleanly if `IBM_QUANTUM_TOKEN` is unset, so
Binder users never hit the live path by accident.

After the run completes the script prints a PASS/FAIL against
[`expected_bounds.json`](expected_bounds.json) (F(K\*) ≥ 0.80,
ΔF ≥ +0.10, ≥ 2 of the runs with ΔF > 0) and exits 0 on PASS, 1 on
FAIL.  The raw JSON lands under `hardware/results/`.

## Offline reanalysis

If the run JSON already exists (or you want to recheck the shipped
reference against the bounds), no QPU time is needed:

```bash
# Recheck your own result:
python hardware/analyze_results.py hardware/results/w_repeat_results_20260421_143011.json

# Recheck the shipped reference:
python hardware/analyze_results.py --reference

# Aggregate across every run under a directory:
python hardware/analyze_results.py --all hardware/results/
```

The analyzer prints the per-run table, mean ± std, a 95 % bootstrap
CI on ΔF, and the same PASS/FAIL bound check.

## Budget & wall time

| parameter | default | cost impact |
| --- | --- | --- |
| circuits per run | 274 (137 K\* + 137 random) | dominant |
| shots per circuit | 1000 | linear |
| runs | 4 | linear |
| total QPU time on a healthy Eagle/Heron backend | ~2 min | well under the 10 min free tier |
| wall-clock time including queue | typically 5–30 min | backend-dependent |

`--shots` and `--n-runs` scale the budget linearly.  On a busy
backend the queue wait dwarfs the QPU time; the script prints every
job ID so you can watch progress on the IBM dashboard.

## Troubleshooting

- **"env var IBM_QUANTUM_TOKEN is unset"** — export the token first
  (see Prerequisites) or pass `--token-env MYVAR`.
- **"No operational backend with enough qubits is available"** —
  your account has no visible backend with ≥ 4 qubits.  Pass
  `--backend ibm_brisbane` (or whatever appears on your dashboard)
  explicitly.
- **Different absolute numbers from the reference** — expected for
  any backend other than `ibm_fez`, or on `ibm_fez` weeks after the
  original run (calibration drift).  The `expected_bounds.json`
  thresholds are set a few σ below the reference means; a PASS there
  is evidence the K\* advantage reproduces, not that the run matches
  to four decimals.
- **Connection dropped mid-run** — a
  `w_repeat_partial_<TIMESTAMP>_run<NN>.json` checkpoint is written
  after each completed run, and every job ID is printed to stdout.
  There is no `--resume` mode: seeds always start at 100, so
  `--n-runs N` reruns seeds `100..`, it does not continue from where
  a prior run stopped.  Two recovery options:
  1. If the partial JSON already has enough completed runs for your
     needs, analyze it directly (pass the specific file — `--all`
     aggregates finalized results, not partials):
     ```
     python hardware/analyze_results.py \
         hardware/results/w_repeat_partial_20260421_143011_run03.json
     ```
  2. Otherwise rerun from scratch with the desired `--n-runs`.  The
     partial JSON and any job IDs in it still serve as provenance.

## Files

| file | purpose |
| --- | --- |
| `run_w_state.py` | main runner: submits circuits, saves JSON |
| `analyze_results.py` | offline analyzer + bound check |
| `circuits.py` | W-state preparation + Pauli-basis rotations |
| `expected_bounds.json` | acceptance thresholds |
| `requirements.txt` | extra pip deps (qiskit, qiskit-ibm-runtime, qiskit-aer) |
| `tests/test_smoke_simulator.py` | local-simulator end-to-end smoke test (no QPU) |
| `tests/test_output_schema.py` | checks emitted JSON layout matches the shipped reference |

## Local-simulator smoke test

Before burning QPU credits you can validate that every code path
(circuit builder, transpile, `SamplerV2(mode=backend)`, counts
extraction, MLE reconstruction, fidelity) works against a local
simulator via `qiskit_ibm_runtime.fake_provider.FakeBrisbane`:

```
pip install -r hardware/requirements.txt
python hardware/tests/test_smoke_simulator.py
```

~6 s wall on a laptop.  Expected output: F(K\*) in the 0.80–0.90
range, ΔF around +0.40 (varies because FakeBrisbane re-samples its
device noise each invocation; not intended to match the flagship
ibm\_fez reference to four decimals).  The test fails if ΔF ≤ 0 —
that invariant only holds when the full pipeline is correct.

The same smoke test runs in CI via
[`../.github/workflows/hardware-smoke.yml`](../.github/workflows/hardware-smoke.yml)
on every push that touches `hardware/`, so qiskit API drift is caught
before it breaks a live-QPU user.  The live-QPU runner itself is
intentionally **not** wired into [`../run_all.py`](../run_all.py) — the
eight built-in tiers replay the shipped hardware snapshot
deterministically, whereas this module spends real QPU credits.
