"""SpatRaster.subset — 0-based indexing (Python convention)."""
import numpy as np
from geospat.rast import rast
from geospat.values import set_values
from geospat._terra import SpatOptions


def test_subset_int_zero_based():
    """Single int is 0-based layer index."""
    r = rast(nrows=2, ncols=2, nlyrs=3)
    r = set_values(r, list(range(2 * 2 * 3)))
    s0 = r.subset(0)
    assert s0.nlyr() == 1
    s1 = r.subset(1)
    assert s1.nlyr() == 1
    assert s0 != s1


def test_subset_list_zero_based():
    """Layer lists use 0-based indices."""
    r = rast(nrows=2, ncols=2, nlyrs=4)
    r = set_values(r, list(range(2 * 2 * 4)))
    a = r.subset([0, 2], SpatOptions())
    assert a.nlyr() == 2


def test_subset_numpy_scalar_int():
    r = rast(nrows=2, ncols=2, nlyrs=2)
    r = set_values(r, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    s = r.subset(np.int64(1))
    assert s.nlyr() == 1


def test_subset_string_layer_name():
    """``subset('name')`` selects the layer whose name matches."""
    r = rast(nrows=2, ncols=2, nlyrs=3)
    r = set_values(r, list(range(2 * 2 * 3)))
    r.names = ["a", "b", "c"]
    sb = r.subset("b")
    assert sb.nlyr() == 1
    assert sb.names == ["b"]
    sa = r.subset("a")
    assert sa.names == ["a"]


def test_subset_negative_index_last_layer():
    """-1 selects last layer (Python-style)."""
    r = rast(nrows=2, ncols=2, nlyrs=3)
    r = set_values(r, list(range(2 * 2 * 3)))
    s = r.subset(-1)
    assert s.nlyr() == 1


def test_getitem_double_bracket_layers():
    """r[[k]] selects layers (0-based)."""
    r = rast(nrows=2, ncols=2, nlyrs=4)
    r = set_values(r, list(range(2 * 2 * 4)))
    a = r[[0]]
    assert a.nlyr() == 1
    b = r[[0, 2]]
    assert b.nlyr() == 2


def test_getitem_single_bracket_cell():
    """r[i] is linear cell index (0-based); returns one value per layer."""
    r = rast(nrows=2, ncols=2, nlyrs=2)
    r = set_values(r, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    v = r[0]
    assert hasattr(v, "shape")
    assert v.shape == (2,)
    assert v[0] == 1.0


def test_getitem_row_col_cell():
    r = rast(nrows=2, ncols=2, nlyrs=1)
    r = set_values(r, [1.0, 2.0, 3.0, 4.0])
    # row 1, col 0 → second row, first column → cell index 2 → value 3
    v = r[1, 0]
    assert float(v[0]) == 3.0


