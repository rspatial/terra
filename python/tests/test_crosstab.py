"""
Tests ported from inst/tinytest/test_crosstab.R

Skipped if crosstab is not yet implemented.
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values

try:
    from terra.generics import crosstab
    _HAS_CROSSTAB = True
except ImportError:
    _HAS_CROSSTAB = False


@pytest.mark.skipif(not _HAS_CROSSTAB, reason="crosstab not yet implemented")
class TestCrosstab:
    def setup_method(self):
        np.random.seed(2)
        r = rast(nrows=10, ncols=10)
        v = np.random.choice(np.arange(1, 6), size=100)
        r = set_values(r, v.astype(float))
        np.random.seed(3)
        s = rast(nrows=10, ncols=10)
        s = set_values(s, np.random.choice(np.arange(1, 4), size=100).astype(float))
        opt = pt.SpatOptions()
        rs = r.deepcopy()
        rs.addSource(s, True, opt)
        self.r = r
        self.s = s
        self.rs = rs

    def test_crosstab_sum_equals_ncell(self):
        ct = crosstab(self.rs)
        assert sum(ct.values.ravel()) == 100

    def test_crosstab_long_ncol(self):
        ct_long = crosstab(self.rs, long=True)
        import pandas as pd
        assert isinstance(ct_long, pd.DataFrame)
        assert ct_long.shape[1] == 3

    def test_crosstab_long_sum(self):
        ct_long = crosstab(self.rs, long=True)
        assert ct_long.iloc[:, 2].sum() == 100
