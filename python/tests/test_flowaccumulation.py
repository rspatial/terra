"""
Tests ported from inst/tinytest/test_flowAccumulation.R
"""
import numpy as np

from geospat.rast import rast
from geospat.extent import ext
from geospat.values import set_values
from geospat.generics import flow_accumulation, terrain


def test_flow_accumulation():
    """Flow direction and accumulation results match expected values."""
    # R: matrix |> rast() uses extent c(0, ncol, 0, nrow), not a global grid.
    elev_vals = [
        78, 72, 69, 71, 58, 49,
        74, 67, 56, 49, 46, 50,
        69, 53, 44, 37, 38, 48,
        64, 58, 55, 22, 31, 24,
        68, 61, 47, 21, 16, 19,
        74, 53, 34, 12, 11, 12,
    ]
    elev = rast(nrows=6, ncols=6, extent=ext(0, 6, 0, 6))
    elev = set_values(elev, np.array(elev_vals, dtype=float))

    flowdir_expected = [
        2, 2, 2, 4, 4, 8,
        2, 2, 2, 4, 4, 8,
        1, 1, 2, 4, 8, 4,
        128, 128, 1, 2, 4, 8,
        2, 2, 1, 4, 4, 4,
        1, 1, 1, 1, 0, 16,
    ]
    flowacc_expected = [
        1, 1, 1, 1, 1, 1,
        1, 2, 2, 3, 3, 1,
        1, 4, 8, 6, 5, 1,
        1, 1, 1, 21, 1, 2,
        1, 1, 1, 2, 25, 1,
        1, 3, 5, 8, 36, 2,
    ]

    flowdir1 = terrain(elev, "flowdir")
    flowacc1 = flow_accumulation(flowdir1)

    fd_vals = np.array(flowdir1.readValues(0, flowdir1.nrow(), 0, flowdir1.ncol()))
    fa_vals = np.array(flowacc1.readValues(0, flowacc1.nrow(), 0, flowacc1.ncol()))

    np.testing.assert_array_equal(fd_vals, flowdir_expected)
    np.testing.assert_array_equal(fa_vals, flowacc_expected)


