"""
Tests ported from inst/tinytest/test_aggregate.R
"""
import math
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values


def _vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def make_r():
    """4×4 raster with values NA, 2..16."""
    r = rast(nrows=4, ncols=4, xmin=0, xmax=1, ymin=0, ymax=1)
    v = [float("nan")] + list(range(2, 17))
    return set_values(r, v)


def _approx_or_nan(a, e):
    if math.isnan(e):
        assert math.isnan(a), f"Expected NaN but got {a}"
    else:
        assert a == pytest.approx(e), f"Expected {e} but got {a}"


def test_aggregate_mean_with_na():
    """aggregate(r, 2, mean) propagates NA."""
    r = make_r()
    f = r.aggregate(2, "mean", na_rm=False)
    result = _vals(f)
    expected = [float("nan"), 5.5, 11.5, 13.5]
    assert len(result) == 4
    for a, e in zip(result, expected):
        _approx_or_nan(a, e)


def test_aggregate_mean_na_rm():
    """aggregate(r, 2, mean, na_rm=True) fills the missing group."""
    r = make_r()
    f = r.aggregate(2, "mean", na_rm=True)
    result = _vals(f)
    np.testing.assert_array_almost_equal(result, [4 + 1/3, 5.5, 11.5, 13.5])


def test_aggregate_multilayer_mean():
    """Two-layer raster: mean without na_rm propagates NA in both layers."""
    r = make_r()
    r2 = r * 2
    opt = pt.SpatOptions()
    rr = r.deepcopy()
    rr.addSource(r2, True, opt)

    f = rr.aggregate(2, "mean", na_rm=False)
    result = _vals(f)
    expected = [float("nan"), 5.5, 11.5, 13.5,
                float("nan"), 11.0, 23.0, 27.0]
    assert len(result) == 8
    for a, e in zip(result, expected):
        _approx_or_nan(a, e)


def test_aggregate_multilayer_mean_na_rm():
    r = make_r()
    r2 = r * 2
    opt = pt.SpatOptions()
    rr = r.deepcopy()
    rr.addSource(r2, True, opt)

    f = rr.aggregate(2, "mean", na_rm=True)
    result = _vals(f)
    np.testing.assert_array_almost_equal(
        result, [4 + 1/3, 5.5, 11.5, 13.5, 8 + 2/3, 11.0, 23.0, 27.0])


def test_aggregate_multilayer_min_na_rm():
    r = make_r()
    r2 = r * 2
    opt = pt.SpatOptions()
    rr = r.deepcopy()
    rr.addSource(r2, True, opt)

    f = rr.aggregate(2, "min", na_rm=True)
    result = _vals(f)
    np.testing.assert_array_almost_equal(result, [2, 3, 9, 11, 4, 6, 18, 22])
