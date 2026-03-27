"""
Tests ported from inst/tinytest/test_focal.R
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra.focal import focal


def _vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def test_focal_asymmetric_weight_matrix():
    """Custom weight matrix with +1 and -1."""
    m = np.zeros((3, 3))
    # R: m[c(4,6)] <- c(1,-1)  (col-major: row=1,col=2 → idx3; row=3,col=2 → idx5)
    m[0, 1] = 1
    m[2, 1] = -1
    r = rast(nrows=3, ncols=3)
    r.setCRS("+proj=merc")
    r = set_values(r, list(range(1, 10)))
    f = focal(r, m, na_rm=True)
    expected = [-4, -5, -6, -6, -6, -6, 4, 5, 6]
    np.testing.assert_array_almost_equal(_vals(f), expected)


def test_focal_transposed_weight_matrix():
    """Transposed weight matrix."""
    m = np.zeros((3, 3))
    m[0, 1] = 1
    m[2, 1] = -1
    r = rast(nrows=3, ncols=3)
    r.setCRS("+proj=merc")
    r = set_values(r, list(range(1, 10)))
    f = focal(r, m.T, na_rm=True)
    expected = [-2, -2, 2, -5, -2, 5, -8, -2, 8]
    np.testing.assert_array_almost_equal(_vals(f), expected)


def test_focal_sum_3x3_ones_na_rm():
    """3×3 ones kernel, na_rm=True."""
    m = np.ones((3, 3))
    r = rast(nrows=3, ncols=3)
    r = set_values(r, list(range(1, 10)))
    f = focal(r, m, na_rm=True)
    expected = [12, 21, 16, 27, 45, 33, 24, 39, 28]
    np.testing.assert_array_almost_equal(_vals(f), expected)


def test_focal_sum_integer_window():
    """Integer window size 3 is equivalent to 3×3 ones kernel."""
    r = rast(nrows=3, ncols=3)
    r = set_values(r, list(range(1, 10)))
    expected = [12, 21, 16, 27, 45, 33, 24, 39, 28]
    f = focal(r, 3, na_rm=True)
    np.testing.assert_array_almost_equal(_vals(f), expected)


def test_focal_no_na_fillvalue():
    """na_rm=False but fillvalue=0 gives same result as na_rm=True for interior."""
    m = np.ones((3, 3))
    r = rast(nrows=3, ncols=3)
    r = set_values(r, list(range(1, 10)))
    f = focal(r, m, na_rm=False, fillvalue=0)
    expected = [12, 21, 16, 27, 45, 33, 24, 39, 28]
    np.testing.assert_array_almost_equal(_vals(f), expected)


def test_focal_no_na_edges_become_nan():
    """na_rm=False with default fillvalue=NaN: only the central cell has a full neighbourhood."""
    m = np.ones((3, 3))
    r = rast(nrows=3, ncols=3)
    r = set_values(r, list(range(1, 10)))
    f = focal(r, m, na_rm=False)
    v = _vals(f)
    assert np.isnan(v[0]) and np.isnan(v[1]) and np.isnan(v[2])  # top row
    assert np.isnan(v[3]) and np.isnan(v[5])                      # left/right middle
    assert np.isnan(v[6]) and np.isnan(v[7]) and np.isnan(v[8])  # bottom row
    assert v[4] == pytest.approx(45)


def test_focal_integer_window_no_na_edges():
    m = np.ones((3, 3))
    r = rast(nrows=3, ncols=3)
    r = set_values(r, list(range(1, 10)))
    # focal(r, 3, na_rm=FALSE) should match focal(r, m, na_rm=FALSE)
    f1 = focal(r, m, na_rm=False)
    f2 = focal(r, 3, na_rm=False)
    v1 = _vals(f1)
    v2 = _vals(f2)
    np.testing.assert_array_equal(
        np.isnan(v1), np.isnan(v2),
        err_msg="NaN pattern should match between matrix and integer-window focal"
    )
    mask = ~np.isnan(v1)
    np.testing.assert_array_almost_equal(v1[mask], v2[mask])
