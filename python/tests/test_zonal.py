"""
Tests ported from inst/tinytest/test_zonal.R

Note: the Python zonal() currently only supports built-in function names
(string arguments). Lambda-based tests are skipped.
"""
import math
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra.zonal import zonal

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False


def make_vz():
    """2×2 rasters: v has a NA, z has zones 1 and 2."""
    v = rast(nrows=2, ncols=2)
    v = set_values(v, [1.0, 2.0, 3.0, float("nan")])
    z = rast(nrows=2, ncols=2)
    z = set_values(z, [1.0, 1.0, 2.0, 2.0])
    return v, z


def _flat(df):
    """
    Flatten a zonal() DataFrame to a list in the same order as
    R's ``unlist(a, use.names=FALSE)`` — zone column first, then value column.
    """
    if isinstance(df, pd.DataFrame):
        return df.values.ravel("F").tolist()
    return list(df)


@pytest.mark.skipif(not _HAS_PANDAS, reason="pandas required")
def test_zonal_isna_string():
    """zonal(v, z, 'isNA') counts NAs per zone."""
    v, z = make_vz()
    a = zonal(v, z, "isNA", na_rm=False)
    flat = _flat(a)
    # zone=1 has cells [1,2] → 0 NAs; zone=2 has [3,NaN] → 1 NA
    # R: c(1,2,0,1)
    assert flat == pytest.approx([1, 2, 0, 1])


@pytest.mark.skipif(not _HAS_PANDAS, reason="pandas required")
def test_zonal_mean_with_na():
    """zonal(v, z, 'mean') with the default na_rm propagates NA."""
    v, z = make_vz()
    a = zonal(v, z, "mean", na_rm=False)
    flat = _flat(a)
    # zone=1 mean([1,2]) = 1.5; zone=2 mean([3,NaN]) = NaN
    assert flat[0] == pytest.approx(1)
    assert flat[1] == pytest.approx(2)
    assert flat[2] == pytest.approx(1.5)
    assert math.isnan(flat[3])


@pytest.mark.skipif(not _HAS_PANDAS, reason="pandas required")
def test_zonal_mean_na_rm_string():
    """zonal(v, z, 'mean', na_rm=True) ignores NAs."""
    v, z = make_vz()
    a = zonal(v, z, "mean", na_rm=True)
    flat = _flat(a)
    assert flat == pytest.approx([1, 2, 1.5, 3])
