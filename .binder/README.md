# Binder configuration

This directory is picked up by [mybinder.org](https://mybinder.org) to
provision an ephemeral Jupyter environment for the K* verification demo.

Files here:
- `runtime.txt` - pins the Python version (3.11).
- `requirements.txt` - pins the four Python packages the notebook needs
  (pyyaml, sympy, numpy, pandas, jupyter).

The Lean toolchain is NOT installed on Binder. The notebook reads
pre-generated artifacts (`lean4/generated/axiom_report.json`,
`verification_manifest.json`) that ship with the repo, so a fresh
Lean build is not required for the 10-minute tier.

For the full 30-minute tier (Lean rebuild), follow the instructions in
[../VERIFICATION.md](../VERIFICATION.md).
