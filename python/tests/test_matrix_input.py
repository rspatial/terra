"""
Tests ported from inst/tinytest/test_matrix-input.R

Creates rasters from numpy arrays using rast() + set_values().
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra._terra import SpatRaster


def make_from_matrix(m, extent=None):
    """
    Create a SpatRaster from a 2-D numpy array, mirroring R's rast(matrix).
    Default extent is (0, ncol, 0, nrow).
    """
    if m.ndim == 2:
        nrows, ncols = m.shape
        nlyrs = 1
        flat = m.ravel(order="C")
    elif m.ndim == 3:
        nrows, ncols, nlyrs = m.shape
        flat = m.reshape(nrows * ncols, nlyrs, order="C").ravel(order="C")
    else:
        raise ValueError("Expected 2-D or 3-D array")

    if extent is None:
        xmin, xmax, ymin, ymax = 0.0, float(ncols), 0.0, float(nrows)
    else:
        xmin, xmax, ymin, ymax = float(extent[0]), float(extent[1]), float(extent[2]), float(extent[3])

    r = rast(nrows=nrows, ncols=ncols, nlyrs=nlyrs,
             xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    r = set_values(r, flat)
    return r


def test_rast_from_2d_array_is_spatraster():
    m = np.arange(1, 13, dtype=float).reshape(3, 4)
    r = make_from_matrix(m)
    assert isinstance(r, SpatRaster)


def test_rast_from_1col_array():
    m = np.arange(1, 11, dtype=float).reshape(10, 1)
    r = make_from_matrix(m)
    assert r.nrow() == 10
    assert r.ncol() == 1
    assert r.nlyr() == 1


def test_rast_from_1row_array():
    m = np.arange(1, 11, dtype=float).reshape(1, 10)
    r = make_from_matrix(m)
    assert r.nrow() == 1
    assert r.ncol() == 10
    assert r.nlyr() == 1


def test_rast_2d_default_extent():
    """Default extent for 3×4 matrix → (0, 4, 0, 3)."""
    m = np.arange(1, 13, dtype=float).reshape(3, 4)
    r = make_from_matrix(m)
    assert r.xmin() == pytest.approx(0)
    assert r.xmax() == pytest.approx(4)
    assert r.ymin() == pytest.approx(0)
    assert r.ymax() == pytest.approx(3)


def test_rast_3d_array_default_extent():
    """3-D array: nrow, ncol from first two dims; nlyr = 3rd dim."""
    a = np.arange(1, 25, dtype=float).reshape(3, 4, 2)
    r = make_from_matrix(a)
    assert r.xmin() == pytest.approx(0)
    assert r.xmax() == pytest.approx(4)
    assert r.ymin() == pytest.approx(0)
    assert r.ymax() == pytest.approx(3)
    assert r.nlyr() == 2


def test_rast_custom_extent():
    """Providing extent overrides the default."""
    m = np.arange(1, 13, dtype=float).reshape(3, 4)
    r = make_from_matrix(m, extent=(-2, 10, -1, 20))
    assert r.xmin() == pytest.approx(-2)
    assert r.xmax() == pytest.approx(10)
    assert r.ymin() == pytest.approx(-1)
    assert r.ymax() == pytest.approx(20)


def test_rast_explicit_dims():
    """Explicit dimension keywords produce correct extent."""
    r = rast(nrows=4, ncols=4, xmin=0, xmax=1, ymin=0, ymax=1)
    assert r.xmin() == pytest.approx(0)
    assert r.xmax() == pytest.approx(1)
    assert r.ymin() == pytest.approx(0)
    assert r.ymax() == pytest.approx(1)
