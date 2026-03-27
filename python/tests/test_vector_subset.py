"""
Tests ported from inst/tinytest/test_vector-subset.R

Requires the bundled lux.shp example file.
"""
import os
import pytest
import numpy as np
import terra as pt
from terra.vect import vect


def find_lux():
    """Return path to lux.shp or skip if not found."""
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "..", "inst", "ex", "lux.shp"),
        r"C:\github\rspatial\terra\inst\ex\lux.shp",
    ]
    for p in candidates:
        p = os.path.normpath(p)
        if os.path.exists(p):
            return p
    pytest.skip("lux.shp not found")


def test_lux_dimensions():
    f = find_lux()
    lux = vect(f)
    assert lux.nrow() == 12
    assert lux.ncol() == 6


def test_lux_id2_values():
    f = find_lux()
    lux = vect(f)
    # R: unlist(lux[,"ID_2",drop=TRUE]) == c(1:7, 12, 8:11)
    expected = list(range(1, 8)) + [12] + list(range(8, 12))
    df = lux.getDF()
    import pandas as pd
    if isinstance(df, dict):
        id2 = list(df["ID_2"])
    else:
        id2 = list(df["ID_2"])
    assert id2 == expected


def test_lux_filter_luxembourg():
    """Filter to Luxembourg rows (NAME_1 == 'Luxembourg')."""
    f = find_lux()
    lux = vect(f)
    df = lux.getDF()
    import pandas as pd
    if isinstance(df, dict):
        df = pd.DataFrame(df)
    idx = list(df.index[df["NAME_1"] == "Luxembourg"])
    x = lux.subset_rows(idx)
    assert x.nrow() == 4
    names = set(pd.DataFrame(x.getDF())["NAME_1"])
    assert names == {"Luxembourg"}
