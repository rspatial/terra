"""
rasterize.py — convert vector data to raster.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages, spatoptions


def _opt() -> SpatOptions:
    return SpatOptions()


_RASTERIZE_FUNS = {"first", "last", "pa", "sum", "mean", "count", "min", "max", "prod"}


def rasterize(
    x: Union[SpatVector, np.ndarray, "pd.DataFrame"],
    y: SpatRaster,
    field: Union[str, float, int, None] = "",
    fun: Optional[Union[str, Callable]] = None,
    *,
    background: float = float("nan"),
    touches: bool = False,
    update: bool = False,
    cover: bool = False,
    by: Optional[str] = None,
    na_rm: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Convert a SpatVector (or coordinate matrix) to a SpatRaster.

    Parameters
    ----------
    x : SpatVector, ndarray (n×2 point coords), or DataFrame
        Features to rasterize.
    y : SpatRaster
        Template raster (extent, resolution, CRS).
    field : str, float, or None
        Attribute column to use as cell values, a numeric constant, or
        ``""`` (presence = 1).
    fun : str or callable, optional
        Aggregation function for overlapping cells.  Built-ins: ``"first"``,
        ``"last"``, ``"sum"``, ``"mean"``, ``"count"``, ``"min"``, ``"max"``,
        ``"pa"`` (presence/absence).
    background : float
        Value for cells not covered by any feature.
    touches : bool
        Include cells touched by polygon boundaries.
    update : bool
        Update (fill) *y* with rasterized values rather than creating a new
        blank raster.
    cover : bool
        Return fractional cell coverage rather than feature values.
    by : str, optional
        Split *x* by the values of this column and rasterize separately
        (returns a multi-layer raster).
    na_rm : bool
        Ignore NA values when aggregating.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    import pandas as pd

    if isinstance(x, np.ndarray) or isinstance(x, pd.DataFrame):
        arr = np.asarray(x, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        x_crds = arr[:, :2]
        values_arr = np.ones(len(x_crds), dtype=float)
        return _rasterize_points_xy(x_crds, values_arr, y, fun, background, update, na_rm, filename, overwrite)

    if not isinstance(x, SpatVector):
        raise TypeError("x must be a SpatVector, ndarray, or DataFrame")

    if by is not None:
        from .generics import split as vect_split
        parts = vect_split(x, by)
        layers = [rasterize(p, y, field=field, fun=fun, background=background,
                            touches=touches, update=update, cover=cover,
                            na_rm=na_rm) for p in parts]
        from .rast import rast as make_rast
        out = make_rast(layers)
        return out

    geom_type = x.geomtype()

    if "points" in geom_type.lower():
        xy = np.array(x.coordinates(), dtype=float).reshape(-1, 2)
        if isinstance(field, (int, float)) and not isinstance(field, bool):
            values_arr = np.full(len(xy), float(field))
        elif isinstance(field, str) and field != "" and field in list(x.names):
            col_idx = list(x.names).index(field)
            vals_raw = [x.getValues1(i, col_idx) for i in range(x.nrow())]
            values_arr = np.array(vals_raw, dtype=float)
        else:
            values_arr = np.ones(len(xy), dtype=float)
        return _rasterize_points_xy(xy, values_arr, y, fun, background, update, na_rm, filename, overwrite)

    # Lines / polygons via C++ rasterize
    if isinstance(field, str) and field != "":
        if field not in list(x.names):
            raise ValueError(f"{field!r} is not a field in x")
    elif field is None or (isinstance(field, str) and field == ""):
        field = ""
    else:
        field = ""

    fun_str = ""
    if fun is not None:
        fun_str = fun if isinstance(fun, str) else getattr(fun, "__name__", "last")

    opt = spatoptions(filename, overwrite)
    xc = y.rasterize(x, field, float(background), touches, fun_str, cover, na_rm, False, cover, opt)
    return messages(xc, "rasterize")


def _rasterize_points_xy(
    xy: np.ndarray,
    values: np.ndarray,
    template: SpatRaster,
    fun: Optional[Union[str, Callable]],
    background: float,
    update: bool,
    na_rm: bool,
    filename: str,
    overwrite: bool,
) -> SpatRaster:
    fun_str = "last"
    if fun is not None:
        fun_str = fun if isinstance(fun, str) else getattr(fun, "__name__", "last")
    if fun_str not in _RASTERIZE_FUNS:
        fun_str = "last"
    opt = spatoptions(filename if not update else "", True if not update else overwrite)
    xc = template.rasterizePointsXY(
        xy[:, 0].tolist(), xy[:, 1].tolist(),
        fun_str, values.tolist(), na_rm, background, opt
    )
    result = messages(xc, "rasterize")
    if update:
        from .generics import cover as rast_cover
        result = rast_cover(result, template, filename=filename, overwrite=overwrite)
    return result


# ---------------------------------------------------------------------------
# rasterizeGeom
# ---------------------------------------------------------------------------

def rasterize_geom(
    x: SpatVector,
    y: SpatRaster,
    fun: str = "count",
    unit: str = "m",
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Rasterize a geometric property (area, length, or count) of *x*.

    Parameters
    ----------
    x : SpatVector
    y : SpatRaster
        Template raster.
    fun : str
        ``"count"`` (default), ``"area"``, or ``"length"``.
    unit : str
        Units for area/length: ``"m"`` or ``"km"``.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    opt = spatoptions(filename, overwrite)
    xc = y.rasterizeGeom(x, unit, fun, opt)
    return messages(xc, "rasterizeGeom")
