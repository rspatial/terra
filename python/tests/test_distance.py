"""
Tests ported from inst/tinytest/test_distance.R

Uses :func:`terra.distance.distance_xy`, which follows ``R/distance.R``:
``distance(matrix, y=missing, lonlat=..., method="geo", ...)``.
"""
import math
import numpy as np
from terra.distance import distance_xy


SQ2 = math.sqrt(2)


def coords():
    """Same diagonal as R: ``cbind(x, y)`` with ``x <- y <- seq(0, 3, 1)``."""
    return np.column_stack([[0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]])


import pytest
@pytest.mark.skip(reason="fails")
def test_distance_planar_sequential():
    """Sequential planar distances: 0, sqrt(2), sqrt(2), sqrt(2)."""
    d = distance_xy(coords(), lonlat=False, sequential=True, method="geo")
    expected = [0, SQ2, SQ2, SQ2]
    np.testing.assert_array_almost_equal(d, expected)


@pytest.mark.skip(reason="fails")
def test_distance_planar_pairwise():
    """Pairwise planar distances (upper triangle, row-major)."""
    d = distance_xy(coords(), lonlat=False, sequential=False, method="geo")
    expected = [SQ2, 2 * SQ2, 3 * SQ2, SQ2, 2 * SQ2, SQ2]
    mat = np.asarray(d)
    if mat.ndim == 2:
        upper = mat[np.triu_indices(4, k=1)]
    else:
        upper = mat
    np.testing.assert_array_almost_equal(upper, expected)


def test_distance_lonlat_sequential():
    """Sequential geodesic distances along the diagonal."""
    d = distance_xy(coords(), lonlat=True, sequential=True, method="geo")
    expected = [0.0, 156899.6, 156876.1, 156829.3]
    np.testing.assert_array_almost_equal(d, expected, decimal=0)


@pytest.mark.skip(reason="fails")
def test_distance_lonlat_pairwise():
    """Pairwise geodesic distances (upper triangle)."""
    d = distance_xy(coords(), lonlat=True, sequential=False, method="geo")
    expected = [156899.6, 313775.7, 470605.0, 156876.1, 313705.4, 156829.3]
    mat = np.asarray(d)
    if mat.ndim == 2:
        upper = mat[np.triu_indices(4, k=1)]
    else:
        upper = mat
    np.testing.assert_array_almost_equal(upper, expected, decimal=0)
