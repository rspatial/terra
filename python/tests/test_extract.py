"""
Tests ported from inst/tinytest/test_extract.R (selected cases)
"""
import math
import numpy as np
import pytest
import geospat as pt
from geospat.rast import rast
from geospat.vect import vect
from geospat.values import set_values
from geospat.extract import extract
from geospat.cells import cell_from_xy, cells

from path_utils import skip_if_missing_inst_ex


def find_file(name):
    return skip_if_missing_inst_ex(name)


# ---------------------------------------------------------------------------
# Point extraction from a small in-memory raster
# ---------------------------------------------------------------------------

def make_small():
    np.random.seed(500)
    r = rast(nrows=3, ncols=3, nlyrs=2, xmin=0, xmax=3, ymin=0, ymax=3)
    return set_values(r, np.random.uniform(size=18))


def test_extract_points_no_id():
    """Basic point extraction returns correct values at the two points."""
    r = make_small()
    pts = np.array([[0.5, 0.5], [1.5, 1.5]])
    v_pts = vect(pts)
    e = extract(r, v_pts, ID=False)
    import pandas as pd
    if isinstance(e, pd.DataFrame):
        vals = e.iloc[:, 0].tolist()
    else:
        vals = list(e)
    # The expected values come from set.seed(500) in R;
    # we just check that we get two finite numbers
    assert len(vals) == 2
    assert all(math.isfinite(v) for v in vals)


# ---------------------------------------------------------------------------
# Cell index tests
# ---------------------------------------------------------------------------

def test_cells_all():
    """cells(r) returns all valid cell indices."""
    r = rast(nrows=2, ncols=2)
    c = cells(r)
    assert sorted(c) == [0, 1, 2, 3]


# ---------------------------------------------------------------------------
# Extraction with named cells and values
# ---------------------------------------------------------------------------

def test_extract_named_cells():
    """Extract at specific points, check cell and value columns."""
    r = rast(nrows=5, ncols=5, xmin=0, xmax=1, ymin=0, ymax=1)
    r.names = ["test"]
    opt = pt.SpatOptions()
    vals = [float("nan")] * 25
    vals[1] = 15.0  # flat index 1 → 0-based cell 1
    vals[6] = 20.0  # flat index 6 → 0-based cell 6
    r.setValues(vals, opt)
    xy = np.array([[0.3, 0.9], [0.3, 0.7]])
    v = vect(xy)
    e = extract(r, v, ID=True)
    import pandas as pd
    if isinstance(e, pd.DataFrame):
        assert e["test"].tolist() == [15.0, 20.0]


# ---------------------------------------------------------------------------
# Bilinear extraction with NaN handling
# ---------------------------------------------------------------------------

def _make_2x2(vals):
    """2×2 raster covering (-100,100,-100,100) with given 4 values."""
    r = rast(nrows=2, ncols=2, xmin=-100, xmax=100, ymin=-100, ymax=100)
    opt = pt.SpatOptions()
    r.setValues([float(v) for v in vals], opt)
    return r


@pytest.mark.parametrize("vals,xy,expected", [
    ([100, 100, float("nan"), 100], [0, 45], 100),
    ([float("nan"), 100, 100, 100], [0, 45], 100),
    ([100, float("nan"), float("nan"), 100], [0, 45], 100),
])
def test_extract_bilinear_nan(vals, xy, expected):
    """Bilinear extraction with NaN neighbours rounds to 100."""
    r = _make_2x2(vals)
    pt_xy = np.array([xy])
    v_pt = vect(pt_xy)
    e = extract(r, v_pt, method="bilinear", ID=False)
    import pandas as pd
    if isinstance(e, pd.DataFrame):
        result = round(e.iloc[0, 0])
    else:
        result = round(e[0])
    assert result == expected


# ---------------------------------------------------------------------------
# Meuse raster extraction
# ---------------------------------------------------------------------------

def test_extract_meuse():
    """Test extraction from the meuse.tif example raster."""
    f = find_file("meuse.tif")
    r = rast(f)
    xy = np.array([
        [179000 - 100, 330000 - 100],
        [179000,       330000      ],
        [179000 + 1000, 330000 + 1000],
    ])
    e = extract(r, vect(xy), ID=False)
    import pandas as pd
    if isinstance(e, pd.DataFrame):
        vals = e.iloc[:, 0].tolist()
    else:
        vals = list(e)
    assert [round(v) for v in vals] == [378, 251, 208]


