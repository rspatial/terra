"""
Tests ported from inst/tinytest/test_merge.R
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra.merge import merge
from terra.names import set_names_inplace


def _make_two_layer(xmin, xmax, val):
    """Create a 2-layer raster covering the given x range."""
    r = rast(nrows=180, ncols=360, xmin=xmin, xmax=xmax, ymin=-90, ymax=90)
    r = set_values(r, float(val))
    opt = pt.SpatOptions()
    r.addSource(r.deepcopy(), True, opt)
    return r


def test_merge_layer_names_preserved():
    """Names from the first raster are preserved in the merge result."""
    r1 = _make_two_layer(0, 1, 1)
    r2 = _make_two_layer(1, 2, 3)
    set_names_inplace(r1, ["x", "y"])
    set_names_inplace(r2, ["x", "y"])

    m = merge(r1, r2)
    assert list(m.names) == ["x", "y"]
