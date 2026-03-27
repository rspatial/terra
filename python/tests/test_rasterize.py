"""
Tests ported from inst/tinytest/test_rasterize.R
"""
import os
import math
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.vect import vect
from terra.rasterize import rasterize


def find_lux():
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "..", "inst", "ex", "lux.shp"),
        r"C:\github\rspatial\terra\inst\ex\lux.shp",
    ]
    for p in candidates:
        p = os.path.normpath(p)
        if os.path.exists(p):
            return p
    pytest.skip("lux.shp not found")


def test_rasterize_cover_by_id():
    """
    Rasterize lux.shp by ID_2 with cover=TRUE.
    Cell (30, 28) (1-based row/col) should have ~1.5% coverage for the
    first polygon and ~98.5% for the fifth.
    """
    f = find_lux()
    v = vect(f)
    r = rast(v, ncols=75, nrows=100)
    z = rasterize(v, r, cover=True, by="ID_2")

    # R: z[30*75+28]  → 1-based cell index
    cell_idx = 30 * 75 + 28 - 1  # convert to 0-based
    nr = z.nrow()
    nc = z.ncol()
    nl = z.nlyr()

    all_vals = np.array(z.readValues(0, nr, 0, nc), dtype=float)
    # Reshape to (ncell, nlyr)
    cell_vals = all_vals.reshape(nc * nr // nc, nc * nl // nl, order='F') \
        if False else all_vals  # keep flat; we need cell_idx-th position per layer

    # With nlyr layers each of ncell cells
    ncell = nr * nc
    vals = [all_vals[layer_i * ncell + cell_idx] for layer_i in range(nl)]

    expected = [0.01538462, float("nan"), float("nan"), float("nan"),
                0.9846154,  float("nan"), float("nan"), float("nan"),
                float("nan"), float("nan"), float("nan"), float("nan")]

    assert len(vals) == len(expected)
    for a, e in zip(vals, expected):
        if math.isnan(e):
            assert math.isnan(a)
        else:
            assert a == pytest.approx(e, rel=2e-5)
