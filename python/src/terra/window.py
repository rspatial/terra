"""
window.py — reading window and raster extension for SpatRaster.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatExtent, SpatOptions
from ._helpers import messages, spatoptions


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# window
# ---------------------------------------------------------------------------

def has_window(x: SpatRaster) -> bool:
    """Return True if *x* has an active reading window."""
    return bool(x.hasWindow())


def set_window(
    x: SpatRaster,
    value: Optional[SpatExtent],
) -> SpatRaster:
    """
    Return a copy of *x* with a reading window applied.

    Parameters
    ----------
    x : SpatRaster
    value : SpatExtent or None
        The window to apply.  Pass None to remove an existing window.

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    if value is None:
        xc.removeWindow()
        return xc

    from .extent import ext as make_ext, intersect_ext
    e = make_ext(x)
    value = intersect_ext(e, value)
    if value is None:
        raise ValueError("window does not overlap with x")
    if not xc.setWindow(value):
        raise RuntimeError("cannot set this window")
    return xc


def remove_window(x: SpatRaster) -> SpatRaster:
    """Return a copy of *x* with any reading window removed."""
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    xc.removeWindow()
    return xc


# ---------------------------------------------------------------------------
# extend (expand raster)
# ---------------------------------------------------------------------------

def extend(
    x: Union[SpatRaster, SpatExtent],
    y: Union[SpatExtent, List[float], float, int],
    *,
    snap: str = "near",
    fill: float = float("nan"),
    filename: str = "",
    overwrite: bool = False,
) -> Union[SpatRaster, SpatExtent]:
    """
    Extend the extent of *x*.

    Parameters
    ----------
    x : SpatRaster or SpatExtent
    y : SpatExtent, or scalar/list
        New extent, or number of cells to add on each side.
        A scalar adds the same amount on all sides.
        A 2-element list [x_cells, y_cells] adds to x and y sides
        respectively.  A 4-element list [left, right, bottom, top] gives
        per-side control.
    snap : str
        Snapping method (``"near"``, ``"in"``, ``"out"``).
    fill : float
        Fill value for added cells.
    filename : str
    overwrite : bool

    Returns
    -------
    Same type as *x*.
    """
    from .extent import ext as make_ext

    if isinstance(x, SpatExtent):
        e_vec = [x.xmin, x.xmax, x.ymin, x.ymax]
        if isinstance(y, (int, float)):
            y = abs(float(y))
            e_vec[0] -= y; e_vec[1] += y
            e_vec[2] -= y; e_vec[3] += y
        elif hasattr(y, '__len__'):
            y_list = [abs(float(v)) for v in y]
            if len(y_list) == 1:
                y_list = y_list * 4
            elif len(y_list) == 2:
                y_list = [y_list[0], y_list[0], y_list[1], y_list[1]]
            elif len(y_list) != 4:
                raise ValueError("y must have 1, 2, or 4 elements")
            e_vec[0] -= y_list[0]; e_vec[1] += y_list[1]
            e_vec[2] -= y_list[2]; e_vec[3] += y_list[3]
        return make_ext(e_vec)

    # SpatRaster
    if not isinstance(y, SpatExtent):
        res_x, res_y = x.xres(), x.yres()
        e_vec = [x.xmin(), x.xmax(), x.ymin(), x.ymax()]
        if isinstance(y, (int, float)):
            y_cells = [abs(int(round(float(y))))] * 4
        else:
            y_cells = [abs(int(round(float(v)))) for v in y]
            if len(y_cells) == 1:
                y_cells = y_cells * 4
            elif len(y_cells) == 2:
                y_cells = [y_cells[0], y_cells[0], y_cells[1], y_cells[1]]
        e_vec[0] -= y_cells[0] * res_x
        e_vec[1] += y_cells[1] * res_x
        e_vec[2] -= y_cells[2] * res_y
        e_vec[3] += y_cells[3] * res_y
        y = make_ext(e_vec)

    opt = spatoptions(filename, overwrite)
    xc = x.expand(y, snap, fill, opt)
    return messages(xc, "extend")
