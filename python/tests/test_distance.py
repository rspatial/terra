"""
Tests ported from inst/tinytest/test_distance.R

R: terra::distance(cbind(x, y), lonlat=FALSE/TRUE, sequential=TRUE/FALSE)
→ Python: distance_vect_self(vect(pts), sequential=..., lonlat=...)
"""
import math
import numpy as np
import pytest
import terra as pt
from terra.vect import vect
from terra.distance import distance_vect_self


SQ2 = math.sqrt(2)


def pts():
    coords = np.column_stack([[0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]])
    return vect(coords)


def test_distance_planar_sequential():
    """Sequential planar distances: 0, sqrt(2), sqrt(2), sqrt(2)."""
    v = pts()
    d = distance_vect_self(v, sequential=True, method="euclidean")
    expected = [0, SQ2, SQ2, SQ2]
    np.testing.assert_array_almost_equal(d, expected)


def test_distance_planar_pairwise():
    """Pairwise planar distances (upper triangle, row-major)."""
    v = pts()
    d = distance_vect_self(v, sequential=False, method="euclidean")
    # Upper triangle (excluding diagonal): (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    expected = [SQ2, 2 * SQ2, 3 * SQ2, SQ2, 2 * SQ2, SQ2]
    # d is returned as a 4×4 matrix; extract upper triangle
    mat = np.asarray(d)
    if mat.ndim == 2:
        upper = mat[np.triu_indices(4, k=1)]
    else:
        upper = mat
    np.testing.assert_array_almost_equal(upper, expected)


def test_distance_lonlat_sequential():
    """Sequential geodesic distances along the diagonal."""
    v = pts()
    d = distance_vect_self(v, sequential=True, method="haversine")
    expected = [0.0, 156899.6, 156876.1, 156829.3]
    np.testing.assert_array_almost_equal(d, expected, decimal=0)


def test_distance_lonlat_pairwise():
    """Pairwise geodesic distances (upper triangle)."""
    v = pts()
    d = distance_vect_self(v, sequential=False, method="haversine")
    expected = [156899.6, 313775.7, 470605.0, 156876.1, 313705.4, 156829.3]
    mat = np.asarray(d)
    if mat.ndim == 2:
        upper = mat[np.triu_indices(4, k=1)]
    else:
        upper = mat
    np.testing.assert_array_almost_equal(upper, expected, decimal=0)
