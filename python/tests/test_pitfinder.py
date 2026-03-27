"""
Tests ported from inst/tinytest/test_pitfinder.R
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.generics import terrain

try:
    from terra.generics import pitfinder
    _HAS_PITFINDER = True
except ImportError:
    _HAS_PITFINDER = False


@pytest.mark.skipif(not _HAS_PITFINDER, reason="pitfinder not yet implemented")
def test_pitfinder():
    """Pit finder detects expected pit cells."""
    dx = dy = 1
    elev = np.zeros((9, 9))
    for r in range(9):
        x = (r - 4) * dx
        for c in range(9):
            y = (c - 4) * dy
            elev[r, c] = 10 + 5 * (x**2 + y**2)

    # double the grid (cbind and rbind)
    elev = np.hstack([elev, elev])
    elev = np.vstack([elev, elev])

    r_elev = rast(elev)
    flowdir = terrain(r_elev, "flowdir")
    pits = pitfinder(flowdir)

    v = np.array(pits.readValues(0, pits.nrow(), 0, pits.ncol()), dtype=float)
    nonzero = np.where(~np.isnan(v) & (v != 0))[0]
    # R 1-based cells 77,86,239,248 → 0-based 76,85,238,247
    expected_cells = {76, 85, 238, 247}
    assert set(nonzero.tolist()) == expected_cells
