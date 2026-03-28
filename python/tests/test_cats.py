"""
Tests ported from inst/tinytest/test_cats.R

Covers: levels / categories on SpatRaster.
"""
import numpy as np
import pandas as pd
import pytest
import geospat as pt
from geospat.rast import rast
from geospat.values import set_values
from geospat.levels import (
    levels, set_levels,
    cats, drop_levels, catalyze,
)


def _vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def make_r():
    np.random.seed(0)
    r = rast(nrows=10, ncols=10)
    v = np.random.choice([1, 2, 3], size=100).astype(float)
    return set_values(r, v)


def test_cats_roundtrip():
    """Categories set with set_levels() can be retrieved with cats()."""
    r = make_r()
    lv = pd.DataFrame({"id": [1, 2, 3, 4],
                       "cover": ["forest", "water", "urban", "other"]})
    r = set_levels(r, lv)
    v = cats(r)[0]
    assert list(v["id"]) == [1, 2, 3, 4]
    assert list(v["cover"]) == ["forest", "water", "urban", "other"]


def test_drop_levels():
    """drop_levels() removes categories not present in the data."""
    r = make_r()
    lv = pd.DataFrame({"id": [1, 2, 3, 4],
                       "cover": ["forest", "water", "urban", "other"]})
    r = set_levels(r, lv)
    r = drop_levels(r)
    remaining = cats(r)[0]
    # values are 1,2,3 only → 'other' (id=4) should be removed
    assert 4 not in list(remaining["id"])
    assert "other" not in list(remaining.iloc[:, 1])


def test_catalyze_nlyr():
    """catalyze() converts a factor raster to numeric layers."""
    np.random.seed(1)
    r = rast(nrows=10, ncols=10)
    r = set_values(r, np.random.choice(np.arange(1, 6), size=100).astype(float))
    lv = pd.DataFrame({"ID": range(1, 6),
                       "val1": range(101, 106),
                       "val2": [v * 0.1 for v in range(1, 6)]})
    r = set_levels(r, lv)
    r_cat = catalyze(r)
    assert r_cat.nlyr() == 2


def test_catalyze_layer_names():
    np.random.seed(1)
    r = rast(nrows=10, ncols=10)
    r = set_values(r, np.random.choice(np.arange(1, 6), size=100).astype(float))
    lv = pd.DataFrame({"ID": range(1, 6),
                       "val1": range(101, 106),
                       "val2": [v * 0.1 for v in range(1, 6)]})
    r = set_levels(r, lv)
    r_cat = catalyze(r)
    names = list(r_cat.names)
    assert "val1" in names
    assert "val2" in names


def test_catalyze_val1_values():
    """val1 layer should equal original integer values + 100."""
    np.random.seed(1)
    r = rast(nrows=10, ncols=10)
    v_orig = np.random.choice(np.arange(1, 6), size=100).astype(float)
    r = set_values(r, v_orig)
    lv = pd.DataFrame({"ID": range(1, 6),
                       "val1": range(101, 106),
                       "val2": [v * 0.1 for v in range(1, 6)]})
    r = set_levels(r, lv)
    r_cat = catalyze(r)
    # Retrieve val1 layer (layer 0)
    opt = pt.SpatOptions()
    r_val1 = r_cat.subset([0], opt)
    v_cat1 = _vals(r_val1)
    # Re-read the original values (may have been modified by set_levels)
    v_r = _vals(r)
    np.testing.assert_array_equal(v_cat1, v_r + 100)


