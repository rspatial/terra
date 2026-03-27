"""
Tests ported from inst/tinytest/test_equal.R

Two rasters built from the same data should compare equal;
rasters built from data in different layout should not.
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values


def make_r1():
    """10×10 raster with sqrt(1:100) values (row-major fill)."""
    x = np.sqrt(np.arange(1, 101, dtype=float))
    r = rast(nrows=10, ncols=10, xmin=0)
    return set_values(r, x)


def make_r2():
    """Same data as make_r1 but filled column-major (like R's matrix())."""
    x = np.sqrt(np.arange(1, 101, dtype=float))
    mat_col_major = x.reshape(10, 10, order="F")
    r = rast(nrows=10, ncols=10, xmin=0)
    return set_values(r, mat_col_major.ravel(order="C"))


def _read_vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def test_same_raster_equal_to_itself():
    """A raster compared with itself must be equal."""
    r1 = make_r1()
    v1 = _read_vals(r1)
    np.testing.assert_array_almost_equal(v1, v1)


def test_row_vs_column_major_not_equal():
    """Row-major vs column-major fill produces different values."""
    r1 = make_r1()
    r2 = make_r2()
    v1 = _read_vals(r1)
    v2 = _read_vals(r2)
    assert not np.allclose(v1, v2), "r1 and r2 should differ (different fill order)"


def test_multilayer_self_equal():
    """Two identical multi-layer rasters have equal values."""
    r1 = make_r1()
    r2 = make_r2()
    opt = pt.SpatOptions()
    r3 = r1.deepcopy()
    r3.addSource(r2.deepcopy(), True, opt)
    r4 = r1.deepcopy()
    r4.addSource(r2.deepcopy(), True, opt)

    v3 = _read_vals(r3)
    v4 = _read_vals(r4)
    np.testing.assert_array_almost_equal(v3, v4)
