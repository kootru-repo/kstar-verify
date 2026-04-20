FROM sagemath/sagemath:latest

WORKDIR /verify
COPY tier1-sage/ /verify/tier1-sage/
COPY tier2-sympy/ /verify/tier2-sympy/
COPY tier3-numpy/ /verify/tier3-numpy/
COPY tier5-fidelity/ /verify/tier5-fidelity/
COPY registry.py /verify/
COPY proofs_registry.yaml /verify/
COPY requirements.txt /verify/

# Install Python dependencies inside SageMath's Python
RUN sage -pip install --no-cache-dir -r requirements.txt 2>/dev/null || \
    pip install --no-cache-dir -r requirements.txt

# Run SageMath tier (44 tests) + Python tiers (2, 3, 5).
# Tiers 4 and 6 use the repo's bundled data/ via validate.sh; run
# `bash validate.sh` from the container for the full suite.
CMD echo "=== Tier 1: SageMath exact arithmetic (44 tests) ===" && \
    cd /verify/tier1-sage && sage run_all.sage --export && \
    echo "" && \
    echo "=== Tier 2: SymPy exact rational arithmetic ===" && \
    cd /verify/tier2-sympy && sage -python run_all.py && \
    echo "" && \
    echo "=== Tier 3: NumPy numerical bounds ===" && \
    cd /verify/tier3-numpy && sage -python run_all.py && \
    echo "" && \
    echo "=== Tier 5: Fidelity HS->F recovery (10 checks) ===" && \
    cd /verify/tier5-fidelity && sage -python test_fidelity_recovery.py && \
    echo "" && \
    echo "Tiers 1, 2, 3, 5 passed (tiers 4, 6 require hardware data — use validate.sh)."
