"""Load and expose claim values from proofs_registry.yaml.

The registry is the single source of truth for all manuscript claims.
Verification scripts compute values independently and compare against
the registry. If a claim changes in proofs.tex, updating the registry
propagates to every test automatically.

Usage:
    from registry import claims
    M = claims.get("thm:basin", "n4_M")          # -> 137
    eps = claims.get("prop:purity_main", "eps_pos_W_k2")  # -> Rational(11, 256)
"""

import os
import ast
import json
from pathlib import Path
from fractions import Fraction


def _parse_value(raw):
    """Convert a registry value string to a Python numeric type.

    Handles: int, float, Fraction (from "11/256"), list, string comparisons.
    """
    if isinstance(raw, (int, float)):
        return raw
    if isinstance(raw, list):
        return [_parse_value(v) for v in raw]

    s = str(raw).strip()

    # Strip trailing comments like "= 0.0430"
    if "=" in s and "/" in s.split("=")[0]:
        s = s.split("=")[0].strip()

    # "~4.9" -> approximate float; "~1/3 for q..." -> string
    if s.startswith("~"):
        rest = s[1:].strip()
        try:
            return float(rest)
        except ValueError:
            return s

    # "> 0.999" -> keep as string (comparison, not value)
    if s.startswith(">") or s.startswith("<"):
        return s

    # "q^2 (empirical)" -> keep as string
    if any(c.isalpha() for c in s) and "/" not in s:
        return s

    # "[9, 56, 24, 32, 16]" -> list of ints
    if s.startswith("["):
        return ast.literal_eval(s)

    # "256/64 = 4" -> parse the fraction part
    if "=" in s:
        s = s.split("=")[0].strip()

    # "11/256" -> Fraction
    if "/" in s:
        try:
            return Fraction(s)
        except (ValueError, ZeroDivisionError):
            return s

    # Plain number
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


class Registry:
    """Loaded claim registry with lookup by statement ID and key."""

    def __init__(self, path):
        # Use a simple YAML parser (avoid pyyaml dependency for portability)
        self._data = self._load_yaml(path)
        self._claims = {}
        for stmt in self._data.get("statements", []):
            sid = stmt["id"]
            kv = stmt.get("key_values", {})
            self._claims[sid] = {k: _parse_value(v) for k, v in kv.items()}
        # Legacy-ID aliases: some tier scripts predate the 2026-04 rename
        # of `prop:*_main` to `lem:*_main`.  Mirror the values so both
        # lookups resolve identically; adding or removing either form
        # updates the other automatically.
        _LEGACY_ALIASES = [
            ("prop:purity_main", "lem:purity_main"),
            ("prop:spectral_q_main", "lem:spectral_q_main"),
        ]
        for old, new in _LEGACY_ALIASES:
            if new in self._claims and old not in self._claims:
                self._claims[old] = self._claims[new]
            elif old in self._claims and new not in self._claims:
                self._claims[new] = self._claims[old]
            elif old in self._claims and new in self._claims:
                # Both forms exist: they must agree, otherwise callers get
                # different values depending on which spelling they use.
                if self._claims[old] != self._claims[new]:
                    raise ValueError(
                        f"registry alias mismatch: '{old}' and '{new}' "
                        f"both present but key_values differ; resolve in "
                        f"proofs_registry.yaml"
                    )

    def get(self, statement_id, key):
        """Get a claim value. Raises KeyError if not found."""
        return self._claims[statement_id][key]

    def statement_ids(self):
        """List all statement IDs."""
        return list(self._claims.keys())

    def keys(self, statement_id):
        """List all key_values for a statement."""
        return list(self._claims.get(statement_id, {}).keys())

    @staticmethod
    def _load_yaml(path):
        """Minimal YAML loader — handles the registry structure without pyyaml."""
        try:
            import yaml
            with open(path, encoding="utf-8") as f:
                return yaml.safe_load(f)
        except ImportError:
            pass
        # Fallback: use pyyaml-free JSON if someone converts, or raise
        json_path = Path(path).with_suffix(".json")
        if json_path.exists():
            with open(json_path) as f:
                return json.load(f)
        raise ImportError(
            f"pyyaml not installed and {json_path} not found. "
            f"Install pyyaml: pip install pyyaml"
        )


def _load_auto():
    """Load registry from KSTAR_REGISTRY env var, or auto-discover next to this file.

    Returns None if no registry can be found.  Downstream callers must
    guard against `claims is None` (or call `registry.load(path)`
    explicitly) before using `claims.get(...)`; otherwise the first
    access raises AttributeError with no context.
    """
    path = os.environ.get("KSTAR_REGISTRY")
    if path is not None:
        return Registry(path)
    # Auto-discover proofs_registry.yaml next to registry.py
    default = Path(__file__).resolve().parent / "proofs_registry.yaml"
    if default.exists():
        return Registry(str(default))
    return None


class _MissingRegistry:
    """Placeholder for the module-level `claims` singleton when auto-load
    fails.  Attribute access raises AttributeError with a descriptive
    message so the error propagates naturally while preserving
    `hasattr()` compatibility.

    Why AttributeError (not RuntimeError): `hasattr(obj, name)` is
    implemented as "call getattr, catch AttributeError, return
    True/False".  Raising RuntimeError would convert a quiet feature
    probe into a crash.  Python's pickling, copying, and repr
    machinery also rely on AttributeError for missing attributes.  By
    raising AttributeError (with a helpful message embedded) we both
    keep those mechanisms working AND surface the root cause the first
    time a caller actually tries to use the registry.

    `bool(claims)` evaluates to False, enabling
    `if claims: claims.get(...)` as a defensive guard.
    """

    __slots__ = ()

    def __getattr__(self, name):
        if name.startswith("__") and name.endswith("__"):
            # Dunder probes: quiet miss for repr/pickle/copy.
            raise AttributeError(name)
        raise AttributeError(
            f"registry.claims was not loaded: attribute {name!r} is "
            f"unavailable.  Set KSTAR_REGISTRY to a proofs_registry.yaml "
            f"path, place the file next to registry.py, or call "
            f"`registry.load(path)` explicitly before using claims."
        )

    def __bool__(self):
        return False

    def __repr__(self):
        return "<registry.claims NOT LOADED; see registry._MissingRegistry>"


# Module-level singleton -- auto-loaded from env or default location.
# If auto-load returns None, wrap it so any downstream access raises a
# descriptive error instead of a bare AttributeError.
claims = _load_auto() or _MissingRegistry()


def load(path):
    """Explicitly load a registry file. Updates the module-level singleton."""
    global claims
    claims = Registry(path)
    return claims
