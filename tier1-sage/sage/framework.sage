"""
SageMath Exact Cross-Verification Framework
============================================

Independent CAS verification of exact-arithmetic claims.
Each domain module registers test functions; the runner dispatches all.

Usage:
    sage run_all.sage                        # full suite
    sage run_all.sage --domain combinatorics # single domain
    sage run_all.sage --export               # write certified_values.json
"""

import json
import os
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Resolve base directory robustly
# SageMath preprocesses .sage -> .py in a temp dir, so __file__ is unreliable.
# ---------------------------------------------------------------------------

def _resolve_base():
    """Find the sage-exact-verification root directory."""
    env_base = os.environ.get("SAGE_VERIFY_BASE")
    if env_base:
        return Path(env_base)
    try:
        candidate = Path(__file__).resolve().parent.parent
        if (candidate / "manifest.yaml").exists():
            return candidate
    except Exception:
        pass
    cwd = Path.cwd()
    for p in [cwd, cwd.parent, cwd.parent.parent]:
        if (p / "manifest.yaml").exists():
            return p
    raise FileNotFoundError(
        "Cannot find sage-exact-verification root (no manifest.yaml). "
        "Run from the sage-exact-verification directory or set SAGE_VERIFY_BASE."
    )

_BASE = _resolve_base()
_SAGE_DIR = _BASE / "sage"

# ---------------------------------------------------------------------------
# Test registry
# ---------------------------------------------------------------------------

_TEST_REGISTRY = {}   # {domain_id: {test_name: callable}}
_CERTIFIED = {}       # {test_id: {value, type, timestamp}}


def register_test(domain, name):
    """Decorator: register a verification test."""
    def decorator(fn):
        if domain not in _TEST_REGISTRY:
            _TEST_REGISTRY[domain] = {}
        if name in _TEST_REGISTRY[domain]:
            raise ValueError(f"Duplicate test registration: {domain}/{name}")
        _TEST_REGISTRY[domain][name] = fn
        return fn
    return decorator


def certify(test_id, value, value_type="exact"):
    """Record a certified value for export. Warns on duplicate test_id."""
    if test_id in _CERTIFIED:
        print(f"  WARNING: certify() overwriting existing value for '{test_id}'")
    _CERTIFIED[test_id] = {
        "value": str(value),
        "type": value_type,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


# ---------------------------------------------------------------------------
# Domain loading
# ---------------------------------------------------------------------------

ALL_DOMAINS = ["combinatorics", "linear_algebra", "probability", "quantum_info", "number_theory"]


def load_domains(domain_filter=None):
    """Load domain modules. If domain_filter given, load only that domain."""
    domains_dir = _SAGE_DIR / "domains"
    targets = [domain_filter] if domain_filter else ALL_DOMAINS
    for d in targets:
        sf = domains_dir / f"{d}.sage"
        if sf.exists():
            load(str(sf))
        else:
            print(f"  WARNING: domain file not found: {sf}")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_test(domain, test_name):
    """Run a single registered test. Returns (passed: bool, message: str)."""
    if domain not in _TEST_REGISTRY:
        return False, f"Domain '{domain}' not loaded"
    if test_name not in _TEST_REGISTRY[domain]:
        return False, f"Test '{test_name}' not in domain '{domain}'"
    fn = _TEST_REGISTRY[domain][test_name]
    try:
        result = fn()
        if result is None:
            return True, "PASS"
        return True, str(result)
    except AssertionError as e:
        return False, f"FAIL: {e}"
    except Exception as e:
        return False, f"ERROR: {type(e).__name__}: {e}"


def run_all(domain_filter=None, export=False):
    """Run all verification tests."""
    load_domains(domain_filter)

    total_tests = sum(len(t) for t in _TEST_REGISTRY.values())
    print(f"\n  SageMath Exact Cross-Verification Suite")
    print(f"  Domains loaded: {sorted(_TEST_REGISTRY.keys())}")
    print(f"  Tests registered: {total_tests}")

    total_pass = 0
    total_fail = 0

    for domain in sorted(_TEST_REGISTRY.keys()):
        tests = _TEST_REGISTRY[domain]
        print(f"\n{'='*60}")
        print(f"  [{domain}] ({len(tests)} tests)")
        print(f"{'='*60}")

        for name in sorted(tests.keys()):
            passed, msg = run_test(domain, name)
            status = "PASS" if passed else "FAIL"
            print(f"  [{status}] {name}")
            if not passed:
                print(f"         {msg}")
            if passed:
                total_pass += 1
            else:
                total_fail += 1

    print(f"\n{'='*60}")
    print(f"  TOTAL: {total_pass} PASS, {total_fail} FAIL, {len(_CERTIFIED)} certified values")
    print(f"{'='*60}")

    if export and _CERTIFIED:
        export_path = _BASE / "certified_values.json"
        with open(export_path, "w") as f:
            json.dump(_CERTIFIED, f, indent=2, sort_keys=True)
        print(f"\n  Exported to {export_path}")

    return total_fail == 0
