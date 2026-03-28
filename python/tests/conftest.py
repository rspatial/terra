"""
Shared pytest fixtures and configuration for the geospat Python test suite.
"""
import os
import sys

import pytest

# This directory must be on sys.path so ``import path_utils`` resolves (pytest
# rootdir is ``python/``, not ``python/tests/``).
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

# Ensure the geospat package is importable when running tests directly
_SRC = os.path.normpath(os.path.join(_TESTS_DIR, "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from path_utils import inst_ex_file


# ---------------------------------------------------------------------------
# Shared helper
# ---------------------------------------------------------------------------

def _find_ex(name):
    """Return absolute path to an example file, or None if not found."""
    p = inst_ex_file(name)
    return p if os.path.exists(p) else None


def skip_if_missing(name):
    p = _find_ex(name)
    if p is None:
        pytest.skip(f"{name} not found in inst/ex/")
    return p
