"""
pitfinder — locate pits (flow sinks) from a flow-direction raster (R ``pitfinder``).
"""
from __future__ import annotations

from ._terra import SpatRaster
from ._helpers import messages, spatoptions


def pitfinder(
    x: SpatRaster,
    *,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Find pit cells from a flow-direction raster.

    Parameters
    ----------
    x : SpatRaster
        Flow directions (e.g. from ``terrain(..., 'flowdir')``).

    Returns
    -------
    SpatRaster
        Non-zero values mark pits (see R ``pitfinder``).
    """
    opt = spatoptions(filename, overwrite)
    xc = x.pitfinder2(opt)
    return messages(xc, "pitfinder")
