"""
Tests ported from inst/tinytest/test_arith.R
"""
import numpy as np
import geospat as pt
from geospat.values import set_values
from geospat.rast import rast


def _vals(r):
    """Return flat float array of all cell values."""
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def make_r():
    """2×2 raster filled with 0.5."""
    r = rast(nrows=2, ncols=2)
    return set_values(r, 0.5)


def test_scalar_minus_raster():
    """2 - r should give 1.5 everywhere."""
    r = make_r()
    result = _vals(2 - r)
    np.testing.assert_array_almost_equal(result, [1.5, 1.5, 1.5, 1.5])


def test_raster_minus_scalar():
    """r - 2 should give -1.5 everywhere."""
    r = make_r()
    result = _vals(r - 2)
    np.testing.assert_array_almost_equal(result, [-1.5, -1.5, -1.5, -1.5])


def test_scalar_divide_raster():
    """2 / r should give 4 everywhere."""
    r = make_r()
    result = _vals(2 / r)
    np.testing.assert_array_almost_equal(result, [4.0, 4.0, 4.0, 4.0])


def test_raster_divide_scalar():
    """r / 2 should give 0.25 everywhere."""
    r = make_r()
    result = _vals(r / 2)
    np.testing.assert_array_almost_equal(result, [0.25, 0.25, 0.25, 0.25])


