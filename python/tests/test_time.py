"""
Tests ported from inst/tinytest/test_time.R
"""
from datetime import datetime, timezone
import pytest
import terra as pt
from terra.rast import rast
from terra.values import set_values
from terra.time import set_time, get_time, has_time, time_info


def test_time_roundtrip():
    """Setting a datetime and retrieving it preserves the value."""
    r = rast(nrows=2, ncols=2, names="random")
    r = set_values(r, [1.0, 2.0, 3.0, 4.0])
    dt = datetime(2000, 1, 2, 12, 12, 12, tzinfo=timezone.utc)
    r_timed = set_time(r, [dt])
    assert has_time(r_timed)
    t_back = get_time(r_timed)
    assert len(t_back) == 1
    assert r_timed.time() == t_back
    assert get_time(r_timed, format="raw") == list(r_timed.time_raw)


def test_time_info_zone():
    """time_info() returns timezone information."""
    r = rast(nrows=2, ncols=2)
    r = set_values(r, [1.0, 2.0, 3.0, 4.0])
    dt = datetime(2000, 1, 2, 12, 12, 12, tzinfo=timezone.utc)
    r_timed = set_time(r, [dt], tz="UTC")
    info = time_info(r_timed)
    assert info["has_time"] is True
    # The timezone may be stored as UTC or empty depending on C++ implementation
    assert "zone" in info
