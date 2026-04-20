"""
SageMath Exact Cross-Verification — Top-level entry point.

Usage:
    sage run_all.sage                        # full suite
    sage run_all.sage --domain combinatorics # single domain
    sage run_all.sage --export               # write certified_values.json
    sage run_all.sage --list                 # list all tests
"""

load("sage/framework.sage")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="SageMath Exact Cross-Verification")
    parser.add_argument("--domain", default=None, help="Run only this domain")
    parser.add_argument("--export", action="store_true", help="Export certified values JSON")
    parser.add_argument("--list", action="store_true", help="List registered tests")
    args = parser.parse_args()

    if args.list:
        load_domains(args.domain)
        for domain, tests in sorted(_TEST_REGISTRY.items()):
            print(f"\n  [{domain}]")
            for name in sorted(tests):
                print(f"    {name}")
        sys.exit(0)

    success = run_all(args.domain, args.export)
    # Use os._exit to bypass Sage's SystemExit handling, which otherwise
    # intercepts sys.exit(0), prints the value, and re-exits with status 1.
    import os
    os._exit(0 if success else 1)
