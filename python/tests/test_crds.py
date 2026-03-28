"""
Tests ported from inst/tinytest/test_crds.R
"""
import math
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.vect import vect
from terra.cells import cell_from_xy


def test_cell_from_xy_valid_origin():
    """(0, 0) → cell 32580 on the default global raster (0-based; R reports 32581)."""
    r = rast()
    cells = cell_from_xy(r, np.array([[0.0, 0.0]]))
    assert cells[0] == 32580


def test_cell_from_xy_all_nan_returns_invalid():
    """NaN coordinate pairs should return an invalid / NA sentinel."""
    r = rast()
    xy = np.array([[float("nan"), float("nan")], [float("nan"), 0.0]])
    cells = cell_from_xy(r, xy)
    # R returns NA; Python may return NaN, -1, or a sentinel.
    # Valid 0-based cells are in [0, ncell - 1].
    ncell = r.ncell()
    for c in cells:
        assert c < 0 or c >= ncell or math.isnan(float(c)), (
            f"Expected invalid cell but got {c}"
        )


def test_vect_point_from_matrix():
    """vect() from a 2-column coordinate matrix creates a point vector."""
    m = np.array([[0, 0], [1, 0]])
    v = vect(m)
    assert v.nrow() == 2
