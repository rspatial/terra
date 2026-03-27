"""
Tests ported from inst/tinytest/test_patches.R
"""
import math
import numpy as np
import pytest
from terra.rast import rast
from terra.values import set_values
from terra.generics import patches


def _vals(r):
    r.readStart()
    try:
        return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)
    finally:
        r.readStop()


def test_patches_d8_zero_as_na():
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [0.0, 1.0, 1.0, 0.0])
    p = patches(r, directions=8, zero_as_na=True)
    v = _vals(p)
    assert math.isnan(v[0]) and math.isnan(v[3])
    # both 1-cells connect diagonally → same patch
    assert v[1] == v[2]


def test_patches_d4_zero_as_na():
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [0.0, 1.0, 1.0, 0.0])
    p = patches(r, directions=4, zero_as_na=True)
    v = _vals(p)
    assert math.isnan(v[0]) and math.isnan(v[3])
    # two 1-cells are only diagonally adjacent → different patches under d=4
    assert v[1] != v[2]


def test_patches_diagonal_d8():
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [1.0, 0.0, 0.0, 1.0])
    p = patches(r, directions=8, zero_as_na=True)
    v = _vals(p)
    assert math.isnan(v[1]) and math.isnan(v[2])
    # diagonal cells connect under d=8 → same patch
    assert v[0] == v[3]


def test_patches_diagonal_d4():
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [1.0, 0.0, 0.0, 1.0])
    p = patches(r, directions=4, zero_as_na=True)
    v = _vals(p)
    assert math.isnan(v[1]) and math.isnan(v[2])
    # under d=4 diagonal is not a connection → different patches
    assert v[0] != v[3]


def test_patches_d4_zero_not_na():
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [1.0, 0.0, 0.0, 1.0])
    p = patches(r, directions=4, zero_as_na=False)
    v = _vals(p)
    # All cells including zeros get a patch id (no NaN values)
    assert all(not math.isnan(v[i]) for i in range(4))
    # All cells connected via zero-cells → one patch
    assert len(set(v.tolist())) == 1


def test_patches_larger_raster_d4():
    """R test_patches.R: 18×36 with four regions → unique patch IDs NaN, 1,2,3,4."""
    r = rast(nrows=18, ncols=36)
    # R: r[1:2,5:8]<-11; r[7:8,1:6]<-12; r[5:6,22:36]<-13; r[15:16,18:29]<-14
    a = np.full((18, 36), np.nan)
    a[0:2, 4:8] = 11
    a[6:8, 0:6] = 12
    a[4:6, 21:36] = 13
    a[14:16, 17:29] = 14
    r = set_values(r, a.ravel(order="C"))
    p = patches(r, directions=4)
    v = _vals(p)
    finite = v[~np.isnan(v)]
    u = np.unique(finite)
    np.testing.assert_array_equal(np.sort(u), [1, 2, 3, 4])
    assert np.isnan(v).any()


def test_patches_3x4_values_d4():
    """
    patches(r, directions=4, values=True) returns patch values.
    Skipped if 'values' kwarg is not supported.

    R ``rast(m)`` uses ``extent = c(0, ncol(m), 0, nrow(m))`` — not the default
    global grid — so the extent must match or patch IDs differ (global wrapping).
    """
    m = np.array([[1, 1, 2, 2],
                  [1, 2, 2, 1],
                  [3, 3, 1, 1]], dtype=float)
    r = rast(nrows=3, ncols=4, xmin=0, xmax=4, ymin=0, ymax=3)
    r = set_values(r, m.ravel(order="C"))
    try:
        p = patches(r, directions=4, values=True)
    except TypeError:
        pytest.skip("patches() does not support values= parameter yet")
    v = _vals(p)
    np.testing.assert_array_equal(v, [1, 1, 2, 2, 1, 2, 2, 4, 5, 5, 4, 4])


def test_patches_3x4_values_d8():
    m = np.array([[1, 1, 2, 2],
                  [1, 2, 2, 1],
                  [3, 3, 1, 1]], dtype=float)
    r = rast(nrows=3, ncols=4, xmin=0, xmax=4, ymin=0, ymax=3)
    r = set_values(r, m.ravel(order="C"))
    try:
        p = patches(r, directions=8, values=True)
    except TypeError:
        pytest.skip("patches() does not support values= parameter yet")
    v = _vals(p)
    np.testing.assert_array_equal(v, [1, 1, 2, 2, 1, 2, 2, 3, 4, 4, 3, 3])
