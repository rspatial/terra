"""
Resolve paths to bundled ``inst/ex/`` example files relative to the repository root.

Tests live in ``python/tests/``; data lives in ``<repo>/inst/ex/``.
"""
from __future__ import annotations

import os

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.normpath(os.path.join(_TESTS_DIR, "..", ".."))


def inst_ex_file(name: str) -> str:
    """Absolute path to ``<repo>/inst/ex/<name>``."""
    return os.path.normpath(os.path.join(_REPO_ROOT, "inst", "ex", name))


def skip_if_missing_inst_ex(name: str) -> str:
    """Return the path to ``inst/ex/name`` or skip the test if the file is absent."""
    import pytest

    p = inst_ex_file(name)
    if not os.path.exists(p):
        pytest.skip(f"{name} not found under inst/ex/")
    return p
