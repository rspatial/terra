"""
Tests ported from inst/tinytest/test_weighted-mean.R
"""
import math
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra.generics import weighted_mean


def _vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def test_weighted_mean_matches_numpy():
    """
    terra::weighted.mean(xt, wt, na.rm=TRUE) should equal
    numpy's manual weighted mean with na_rm.
    """
    # Build single-layer rasters with uniform values
    x = rast(nrows=10, ncols=10, xmin=0, xmax=10)
    x = set_values(x, 1.0)
    y = rast(nrows=10, ncols=10, xmin=0, xmax=10)
    y = set_values(y, 2.0)
    z = rast(nrows=10, ncols=10, xmin=0, xmax=10)
    z = set_values(z, 3.0)

    # weights raster: layers with values 1, 2, 3
    opt = pt.SpatOptions()
    wt = x.deepcopy()
    wt.addSource(y.deepcopy(), True, opt)
    wt.addSource(z.deepcopy(), True, opt)

    # data raster: layers 1, 2, NA-at-cell-1
    xt = x.deepcopy()
    xt.addSource(y.deepcopy(), True, opt)
    z_na = z.deepcopy()
    vals = _vals(z_na).copy()
    vals[0] = float("nan")
    z_na.setValues(vals.tolist(), opt)
    xt.addSource(z_na, True, opt)

    wm = weighted_mean(xt, wt, na_rm=True)
    terra_wm = float(np.nanmin(_vals(wm)))

    # stats::weighted.mean(x=c(1,2,NA), w=1:3, na.rm=TRUE) = (1*1 + 2*2)/(1+2) = 5/3
    stats_wm = (1 * 1 + 2 * 2) / (1 + 2)

    assert terra_wm == pytest.approx(stats_wm, rel=1e-5)
