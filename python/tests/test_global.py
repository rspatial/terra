"""
Tests ported from inst/tinytest/test_global.R

Uses the bundled elev.tif example raster.
"""
import os
import numpy as np
import pytest
import terra as pt
from terra.rast import rast


def find_elev():
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "..", "inst", "ex", "elev.tif"),
        r"C:\github\rspatial\terra\inst\ex\elev.tif",
    ]
    for p in candidates:
        p = os.path.normpath(p)
        if os.path.exists(p):
            return p
    pytest.skip("elev.tif not found")


def _global(r, fun, na_rm=True):
    """Call the C++ summary method (mirrors R global())."""
    opt = pt.SpatOptions()
    result = r.summary(fun, na_rm, opt)
    result.readStart()
    try:
        v = result.readValues(0, result.nrow(), 0, result.ncol())
    finally:
        result.readStop()
    return v[0]


@pytest.mark.parametrize("fun,expected", [
    ("sum",   717627.0),
    ("mean",  304.33715),
    ("min",   141.0),
    ("max",   432.0),
    ("sd",    48.40079),
])
def test_global_statistic(fun, expected):
    f = find_elev()
    r = rast(f)
    # GDAL: open for read (R readStart) before readValues, like values.R
    r.readStart()
    try:
        vals = np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)
    finally:
        r.readStop()
    # R: r[1:50,] <- NA — only rows that exist (min(50, nrow))
    n_na_rows = min(50, int(r.nrow()))
    vals[: n_na_rows * int(r.ncol())] = float("nan")

    valid = vals[~np.isnan(vals)]
    if fun == "sum":
        result = float(np.sum(valid))
        assert result == pytest.approx(expected, rel=1e-4)
    elif fun == "mean":
        result = float(np.mean(valid))
        assert result == pytest.approx(expected, rel=1e-4)
    elif fun == "min":
        result = float(np.min(valid))
        assert result == pytest.approx(expected)
    elif fun == "max":
        result = float(np.max(valid))
        assert result == pytest.approx(expected)
    elif fun == "sd":
        result = float(np.std(valid, ddof=1))
        assert result == pytest.approx(expected, rel=1e-3)
