"""
Tests ported from inst/tinytest/test_geom.R

Covers: resolution, dimensions, extent, and cell-to-xy conversion.
"""
import math
import numpy as np
import pytest
import geospat as pt
from geospat.rast import rast
from geospat.values import set_values
from geospat.cells import xy_from_cell


def make_r3():
    """
    Three-layer 10×20 raster with values 200..1 in the first layer,
    and r/2 and r*2 as the second and third layers.
    """
    r = rast(nrows=10, ncols=20, xmin=-10, xmax=10, ymin=-5, ymax=6)
    r = set_values(r, list(range(200, 0, -1)))
    r2 = r * 2
    r_half = r / 2      # Python operator via arith.py
    # Note: r2 is r*2, r_half is r/2 — but we need r, r/2, r*2 as layers
    opt = pt.SpatOptions()
    out = r.deepcopy()
    out.addSource(r_half, True, opt)
    out.addSource(r2, True, opt)
    return out


def test_resolution():
    r = make_r3()
    assert r.xres() == pytest.approx(1.0)
    assert r.yres() == pytest.approx(1.1)


def test_dimensions():
    r = make_r3()
    assert r.nrow() == 10
    assert r.ncol() == 20
    assert r.nlyr() == 3


def test_extent():
    r = make_r3()
    assert r.xmin() == pytest.approx(-10)
    assert r.xmax() == pytest.approx(10)
    assert r.ymin() == pytest.approx(-5)
    assert r.ymax() == pytest.approx(6)


def test_xy_from_cell_10():
    """0-based cell 9 (same location as R / terra cell 10) → x=-0.5, y=5.45."""
    r = make_r3()
    xy = xy_from_cell(r, 9)
    assert xy[0][0] == pytest.approx(-0.5)
    assert xy[0][1] == pytest.approx(5.45)


def test_xy_from_cell_1():
    """First cell (0-based index 0, top-left centre) → x=-9.5, y=5.45."""
    r = make_r3()
    xy = xy_from_cell(r, 0)
    assert xy[0][0] == pytest.approx(-9.5)
    assert xy[0][1] == pytest.approx(5.45)


def test_xy_from_cell_last():
    """Last cell → x=9.5, y=-4.45."""
    r = make_r3()
    ncell = r.ncell()
    xy = xy_from_cell(r, ncell - 1)
    assert xy[0][0] == pytest.approx(9.5)
    assert xy[0][1] == pytest.approx(-4.45, rel=1e-4)


def test_xy_from_cell_out_of_bounds():
    """Cell -1 or ncell → NaN coords (or out-of-range indicator)."""
    r = make_r3()
    ncell = r.ncell()
    for bad_cell in (-1, ncell):
        xy = xy_from_cell(r, bad_cell)
        # Out-of-range 0-based indices; check coords are NaN or extreme.
        x_coord = float(xy[0][0])
        y_coord = float(xy[0][1])
        is_invalid = (
            math.isnan(x_coord) or math.isnan(y_coord)
            or abs(x_coord) > 1e15 or abs(y_coord) > 1e15
        )
        assert is_invalid, (
            f"Expected invalid coords for cell {bad_cell}, got ({x_coord}, {y_coord})"
        )


