"""
Tests ported from inst/tinytest/test_freq.R

Skipped when freq() is not yet in the Python API.
"""
from collections import Counter
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values

try:
    from terra.generics import freq
    _HAS_FREQ = True
except ImportError:
    _HAS_FREQ = False


@pytest.mark.skipif(not _HAS_FREQ, reason="freq() not yet implemented in Python terra")
class TestFreq:
    def setup_method(self):
        np.random.seed(2)
        r = rast(nrows=10, ncols=10)
        vals = np.random.choice(np.arange(1, 6), size=100)
        self.r = set_values(r, vals.astype(float))
        self.vals = vals

    def test_freq_returns_dataframe(self):
        import pandas as pd
        f = freq(self.r)
        assert isinstance(f, pd.DataFrame)

    def test_freq_columns(self):
        f = freq(self.r)
        assert list(f.columns) == ["layer", "value", "count"]

    def test_freq_counts_match(self):
        f = freq(self.r)
        tbl = Counter(self.vals.tolist())
        f_sorted = f.sort_values("value").reset_index(drop=True)
        expected_vals = sorted(tbl.keys())
        expected_counts = [tbl[v] for v in expected_vals]
        np.testing.assert_array_equal(f_sorted["value"].astype(int).tolist(), expected_vals)
        np.testing.assert_array_equal(f_sorted["count"].tolist(), expected_counts)

    def test_freq_value_filter(self):
        tbl = Counter(self.vals.tolist())
        f = freq(self.r, value=5)
        assert f["count"].iloc[0] == tbl[5]
