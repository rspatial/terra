"""
Shared pytest fixtures and configuration for the terra Python test suite.
"""
import sys
import os
import pytest
import numpy as np

# Ensure the terra package is importable when running tests directly
_SRC = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


# ---------------------------------------------------------------------------
# Shared helper
# ---------------------------------------------------------------------------

def _find_ex(name):
    """Return absolute path to an example file, or None if not found."""
    candidates = [
        os.path.normpath(os.path.join(
            os.path.dirname(__file__), "..", "..", "inst", "ex", name
        )),
        os.path.normpath(os.path.join(
            r"C:\github\rspatial\terra\inst\ex", name
        )),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def skip_if_missing(name):
    p = _find_ex(name)
    if p is None:
        pytest.skip(f"{name} not found in inst/ex/")
    return p
