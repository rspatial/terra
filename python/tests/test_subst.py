"""
Tests ported from inst/tinytest/test_subst.R
"""
import math
import numpy as np
import pytest
import geospat as pt
from geospat.rast import rast
from geospat.values import set_values
from geospat.generics import subst


def _vals(r):
    return np.array(r.readValues(0, r.nrow(), 0, r.ncol()), dtype=float)


def test_subst_integer_lookup():
    """Integer lookup: 1→21, 2→22, values 3-5 unchanged."""
    v = np.tile([1., 2., 3., 4., 5.], 20)[:100]
    r = rast(nrows=10, ncols=10)
    r = set_values(r, v)
    r_sub = subst(r, [1, 2], [21, 22])
    v_sub = _vals(r_sub)
    v_in  = _vals(r)
    np.testing.assert_array_equal(v_sub[v_in == 1], np.full(int((v_in == 1).sum()), 21))
    np.testing.assert_array_equal(v_sub[v_in == 2], np.full(int((v_in == 2).sum()), 22))
    np.testing.assert_array_equal(v_sub[v_in > 2], v_in[v_in > 2])


def test_subst_float_lookup():
    """Float lookup: 0.1→99.1, 0.5→99.5, others unchanged."""
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [0.1, 0.2, 0.5, 0.9])
    r_sub = subst(r, [0.1, 0.5], [99.1, 99.5])
    v_sub = _vals(r_sub)
    v_in  = _vals(r)
    assert v_sub[v_in == 0.1][0] == pytest.approx(99.1)
    assert v_sub[v_in == 0.5][0] == pytest.approx(99.5)
    assert v_sub[v_in == 0.2][0] == pytest.approx(0.2)
    assert v_sub[v_in == 0.9][0] == pytest.approx(0.9)


def test_subst_others():
    """Values not in lookup → replaced with 'others' value (1.5)."""
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [1.0, 2.0, 3.0, 4.0])
    r_sub = subst(r, [1, 2], [11, 12], others=1.5)
    v_sub = _vals(r_sub)
    np.testing.assert_array_almost_equal(v_sub, [11, 12, 1.5, 1.5])


def test_subst_to_na():
    """Substitute values 2-10 with NA, leave 1 unchanged."""
    r = rast(nrows=1, ncols=10)
    r = set_values(r, list(range(1, 11)))
    r_sub = subst(r, list(range(2, 11)), [float("nan")] * 9)
    v_sub = _vals(r_sub)
    assert v_sub[0] == pytest.approx(1)
    assert all(math.isnan(v) for v in v_sub[1:])


