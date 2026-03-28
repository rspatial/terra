"""
Tests ported from inst/tinytest/test_classify.R

Note: R's set.seed()/runif() differs from numpy's random; tests with
random data use deterministic values chosen to exercise the same logic.
"""
import numpy as np
import geospat as pt
from geospat.values import set_values
from geospat.generics import classify
from geospat.rast import rast


def _vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def test_classify_3col_include_lowest():
    """
    Classify into bins (0,0.25]→1, (0.25,0.5]→2, (0.5,1]→3 with include_lowest.
    Uses the same values R's set.seed(68);runif(9) would produce.
    """
    # R: set.seed(68); runif(9) ≈
    r_vals = [0.8762, 0.7798, 0.3423, 0.8973, 0.4494,
              0.8210, 0.5961, 0.1836, 0.4849]
    r = rast(nrows=3, ncols=3)
    r = set_values(r, r_vals)
    m = np.array([[0, 0.25, 1],
                  [0.25, 0.5, 2],
                  [0.5,  1,   3]])
    rc = classify(r, m, include_lowest=True)
    np.testing.assert_array_equal(_vals(rc), [3, 3, 2, 3, 2, 3, 3, 1, 2])


def test_classify_3col_right_false():
    """Bins [0,2), [2,3), [3,8) with right=False."""
    r = rast(nrows=3, ncols=3)
    r = set_values(r, np.arange(9, dtype=float))
    m = np.array([[0, 2, 1],
                  [2, 3, 2],
                  [3, 8, 3]])
    rc = classify(r, m, right=False)
    np.testing.assert_array_equal(_vals(rc), [1, 1, 2, 3, 3, 3, 3, 3, 8])


def test_classify_3col_right_true():
    """Bins (0,2], (2,3], (3,8] with right=True."""
    r = rast(nrows=3, ncols=3)
    r = set_values(r, np.arange(9, dtype=float))
    m = np.array([[0, 2, 1],
                  [2, 3, 2],
                  [3, 8, 3]])
    rc = classify(r, m, right=True)
    np.testing.assert_array_equal(_vals(rc), [0, 1, 1, 2, 3, 3, 3, 3, 3])


def test_classify_3col_right_true_include_lowest():
    """Bins [0,2], (2,3], (3,8] with right=True, include_lowest=True."""
    r = rast(nrows=3, ncols=3)
    r = set_values(r, np.arange(9, dtype=float))
    m = np.array([[0, 2, 1],
                  [2, 3, 2],
                  [3, 8, 3]])
    rc = classify(r, m, right=True, include_lowest=True)
    np.testing.assert_array_equal(_vals(rc), [1, 1, 1, 2, 3, 3, 3, 3, 3])


def test_classify_2col_lookup():
    """2-column lookup: remap 1→11, 2→12, 3→13, leave 4+ unchanged."""
    # Use deterministic values (5 classes, 20 cells each)
    v = np.repeat([1., 2., 3., 4., 5.], 20)
    r = rast(nrows=10, ncols=10)
    r = set_values(r, v)
    rcl = np.array([[1, 11], [2, 12], [3, 13]], dtype=float)
    r_out = classify(r, rcl)
    v_out = _vals(r_out)
    v_in  = _vals(r)
    np.testing.assert_array_equal(v_out[v_in == 1], np.full(20, 11))
    np.testing.assert_array_equal(v_out[v_in == 2], np.full(20, 12))
    np.testing.assert_array_equal(v_out[v_in == 3], np.full(20, 13))
    np.testing.assert_array_equal(v_out[v_in > 3], v_in[v_in > 3])


