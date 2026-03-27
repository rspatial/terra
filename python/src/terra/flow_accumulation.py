"""
flow_accumulation — upslope contributing area from flow directions (R ``flowAccumulation``).
"""
from __future__ import annotations

from typing import Optional

from ._terra import SpatRaster
from ._helpers import messages, spatoptions


def flow_accumulation(
    x: SpatRaster,
    weight: Optional[SpatRaster] = None,
    *,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Compute flow accumulation from a flow-direction raster.

    Parameters
    ----------
    x : SpatRaster
        Flow directions (e.g. from ``terrain(x, "flowdir")``).
    weight : SpatRaster, optional
        If provided, weighted flow accumulation (``flowAccu2_weight``).

    Returns
    -------
    SpatRaster
        Flow accumulation (number of upslope cells, or weighted sum).
    """
    opt = spatoptions(filename, overwrite)
    if weight is None:
        xc = x.flowAccu2(opt)
    else:
        xc = x.flowAccu2_weight(weight, opt)
    return messages(xc, "flowAccumulation")
