FROM sagemath/sagemath:latest

WORKDIR /verify
# `sage` runs as a non-root user inside this base image; `sage-preparse`
# writes temp files next to the .sage source, and every Python tier
# may write __pycache__/.  COPY --chown ensures the sage user owns the
# tree so those writes succeed without a separate RUN chown layer.
COPY --chown=sage:sage tier1-sage/ /verify/tier1-sage/
COPY --chown=sage:sage tier2-sympy/ /verify/tier2-sympy/
COPY --chown=sage:sage tier3-numpy/ /verify/tier3-numpy/
COPY --chown=sage:sage tier5-fidelity/ /verify/tier5-fidelity/
# tier4-independent provides the shared core.py / robust_mle.py
# modules that tier 3 (basis-cover ILP) and tier 5 (fidelity recovery)
# import via sys.path injection.  ~360 KB.
COPY --chown=sage:sage tier4-independent/ /verify/tier4-independent/
# Lean sources (not the .lake build cache) are included so the Python
# tiers can SHA-check the formalisation against the registry without
# needing elan/lake installed in this container.  ~450 KB, 38 files.
COPY --chown=sage:sage lean4/KstarFormal/ /verify/lean4/KstarFormal/
# Small helper scripts imported by the tier runners (discover_tests.py
# provides the import-time test-coverage check used by tier 2 and
# tier 5; lint_commit_messages.py ships for hermeticity).
COPY --chown=sage:sage scripts/ /verify/scripts/
COPY --chown=sage:sage registry.py /verify/
COPY --chown=sage:sage proofs_registry.yaml /verify/
COPY --chown=sage:sage requirements.txt /verify/
COPY --chown=sage:sage run_docker_suite.sh /verify/run_docker_suite.sh
RUN chmod +x /verify/run_docker_suite.sh

# Install Python dependencies inside SageMath's Python
RUN sage -pip install --no-cache-dir -r requirements.txt 2>/dev/null || \
    pip install --no-cache-dir -r requirements.txt

# Run SageMath tier (44 tests) + Python tiers (2, 3, 5).
# Tiers 4 and 6 use the repo's bundled data/ via validate.sh; run
# `bash validate.sh` from the container for the full suite.
CMD ["bash", "/verify/run_docker_suite.sh"]
