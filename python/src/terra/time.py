"""
time.py — time metadata for SpatRaster layers.
"""
from __future__ import annotations
from datetime import datetime, date, timezone
from typing import List, Optional, Union

from ._terra import SpatRaster


# ---------------------------------------------------------------------------
# has_time / time_info
# ---------------------------------------------------------------------------

def has_time(x: SpatRaster) -> bool:
    """Return True if *x* has time metadata."""
    return bool(x.hasTime)


def time_info(x: SpatRaster) -> dict:
    """
    Return a dict describing the time metadata of *x*.

    Returns
    -------
    dict with keys ``"has_time"``, ``"step"``, and ``"zone"``.
    """
    has = bool(x.hasTime)
    if has:
        step = str(x.timestep) if hasattr(x, "timestep") else ""
        zone = str(x.timezone) if hasattr(x, "timezone") else "UTC"
        if step != "seconds":
            zone = ""
    else:
        step, zone = "", ""
    return {"has_time": has, "step": step, "zone": zone}


# ---------------------------------------------------------------------------
# get time
# ---------------------------------------------------------------------------

def get_time(
    x: SpatRaster,
    format: str = "",
) -> Union[List[datetime], List[date], List[float], List[int]]:
    """
    Return the time values associated with each layer.

    Parameters
    ----------
    x : SpatRaster
    format : str
        How to format the output:

        - ``""`` — use the stored time step (default).
        - ``"seconds"`` — return datetime objects.
        - ``"days"`` — return date objects.
        - ``"months"`` — return month integers (1–12).
        - ``"years"`` — return year integers.
        - ``"yearmonths"`` — return year + (month-1)/12 floats.
        - ``"raw"`` — return raw numeric seconds-since-epoch.

    Returns
    -------
    list
    """
    if not x.hasTime:
        return [None] * x.nlyr()

    raw = list(x.time)
    tstep = str(x.timestep) if hasattr(x, "timestep") else "seconds"
    tz_str = str(x.timezone) if hasattr(x, "timezone") else "UTC"
    if not tz_str:
        tz_str = "UTC"

    if format == "" or format is None:
        format = tstep
    if format == "raw":
        return list(raw)

    # Convert raw seconds → datetime objects
    def _to_dt(v):
        try:
            return datetime.fromtimestamp(float(v), tz=timezone.utc)
        except (OSError, OverflowError, ValueError):
            return None

    dts = [_to_dt(v) for v in raw]

    if format == "seconds":
        return dts
    elif format == "days":
        return [d.date() if d is not None else None for d in dts]
    elif format == "months":
        return [d.month if d is not None else None for d in dts]
    elif format == "years":
        return [d.year if d is not None else None for d in dts]
    elif format == "yearmonths":
        return [
            d.year + (d.month - 1) / 12.0 if d is not None else None
            for d in dts
        ]
    else:
        return dts


# ---------------------------------------------------------------------------
# set time
# ---------------------------------------------------------------------------

def set_time(
    x: SpatRaster,
    value: Optional[Union[List, None]],
    tstep: str = "",
    tz: str = "UTC",
) -> SpatRaster:
    """
    Set time metadata on a copy of *x*.

    Parameters
    ----------
    x : SpatRaster
    value : list of datetime, date, or numeric, or None
        Time values for each layer.  None removes time metadata.
    tstep : str
        Time step: ``"seconds"``, ``"days"``, ``"months"``, ``"years"``,
        ``"yearmonths"``, or ``""`` (auto).
    tz : str
        Timezone string (e.g. ``"UTC"``).

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x

    if value is None:
        xc.setTime([], "remove", "")
        return xc

    if len(value) != xc.nlyr():
        raise ValueError("len(value) must equal nlyr(x)")

    # Determine tstep and convert values to seconds-since-epoch
    detected_step = tstep
    epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
    seconds = []
    for v in value:
        if isinstance(v, datetime):
            if v.tzinfo is None:
                v = v.replace(tzinfo=timezone.utc)
            if not detected_step:
                detected_step = "seconds"
            seconds.append((v - epoch).total_seconds())
        elif isinstance(v, date):
            dt = datetime(v.year, v.month, v.day, tzinfo=timezone.utc)
            if not detected_step:
                detected_step = "days"
            seconds.append((dt - epoch).total_seconds())
        elif isinstance(v, (int, float)):
            seconds.append(float(v))
        else:
            seconds.append(float("nan"))

    if not detected_step:
        detected_step = "seconds"

    # C++ setTime(std::vector<int64_t>, step, zone) — step must be one of
    # seconds|raw|days|yearmonths|years|months (not "").
    seconds_i: List[int] = []
    for s in seconds:
        if isinstance(s, float) and s != s:
            seconds_i.append(0)
        else:
            seconds_i.append(int(round(float(s))))

    if not xc.setTime(seconds_i, detected_step, tz):
        raise RuntimeError("could not set time metadata")
    return xc
