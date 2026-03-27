"""
Tests ported from inst/tinytest/test_extent.R
"""
import terra as pt


def test_intersect_disjoint_extents_is_none():
    """Intersection of disjoint extents returns an invalid / empty extent."""
    a = pt.ext(0, 10, 0, 10)
    b = pt.ext(100, 101, 100, 101)
    result = a.intersect(b)
    # R returns NULL; the C++ method returns an invalid SpatExtent
    assert not result.valid
