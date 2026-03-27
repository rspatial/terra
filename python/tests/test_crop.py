"""
Tests ported from inst/tinytest/test_crop.R
"""
import numpy as np
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra.generics import crop
from terra.extent import ext


def _ext_vec(r):
    """Return [xmin, xmax, ymin, ymax] from a raster."""
    return [r.xmin(), r.xmax(), r.ymin(), r.ymax()]


def make_r():
    return rast(nrows=10, ncols=10, xmin=0, xmax=10, ymin=0, ymax=10)


def test_crop_larger_extent():
    r = make_r()
    e = ext(-10, 15, -10, 15)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e)), [0, 10, 0, 10])


def test_crop_larger_extent_extend():
    r = make_r()
    e = ext(-10, 15, -10, 15)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e, extend=True)), [-10, 15, -10, 15])


def test_crop_partial_overlap():
    r = make_r()
    e = ext(-10, 5, -10, 5)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e)), [0, 5, 0, 5])


def test_crop_partial_overlap_extend():
    r = make_r()
    e = ext(-10, 5, -10, 5)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e, extend=True)), [-10, 5, -10, 5])


def test_crop_exact_fit():
    r = make_r()
    e = ext(0, 5, 0, 5)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e)), [0, 5, 0, 5])


def test_crop_exact_fit_extend():
    r = make_r()
    e = ext(0, 5, 0, 5)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e, extend=True)), [0, 5, 0, 5])


def test_crop_with_values_larger():
    r = make_r()
    r = set_values(r, 1.0)
    e = ext(-10, 15, -10, 15)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e)), [0, 10, 0, 10])
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e, extend=True)), [-10, 15, -10, 15])


def test_crop_with_values_partial():
    r = make_r()
    r = set_values(r, 1.0)
    e = ext(-10, 5, -10, 5)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e)), [0, 5, 0, 5])
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e, extend=True)), [-10, 5, -10, 5])


def test_crop_with_values_exact():
    r = make_r()
    r = set_values(r, 1.0)
    e = ext(0, 5, 0, 5)
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e)), [0, 5, 0, 5])
    np.testing.assert_array_almost_equal(_ext_vec(crop(r, e, extend=True)), [0, 5, 0, 5])
