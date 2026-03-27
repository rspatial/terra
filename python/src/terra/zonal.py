"""
zonal.py — zonal statistics.
"""
from __future__ import annotations
from typing import Callable, Optional, Union

from ._terra import SpatRaster, SpatOptions
from ._helpers import messages, spatoptions

_cpp_zonal = SpatRaster.zonal  # captured before any monkey-patching


def _opt() -> SpatOptions:
    return SpatOptions()


def _empty_group_raster(x: SpatRaster) -> SpatRaster:
    """
    Third argument to C++ ``zonal(z, g, ...)`` when there is no grouping layer.

    Matches R ``grast <- rast()`` passed as ``x@pntr$zonal(z@pntr, grast@pntr, ...)``.
    """
    from .rast import rast

    e = x.extent
    v = e.vector
    cr = x.get_crs("wkt")
    return rast(
        None,
        nrows=x.nrow(),
        ncols=x.ncol(),
        nlyrs=1,
        xmin=float(v[0]),
        xmax=float(v[1]),
        ymin=float(v[2]),
        ymax=float(v[3]),
        crs=cr if cr else None,
    )


_ZONAL_FUNS = {
    "sum", "mean", "median", "modal", "min", "max", "prod", "any", "all",
    "count", "sd", "std", "first", "isNA", "notNA",
}


def zonal(
    x: SpatRaster,
    z: SpatRaster,
    fun: Union[str, Callable] = "mean",
    *,
    na_rm: bool = True,
    as_raster: bool = False,
    wide: bool = True,
    filename: str = "",
    overwrite: bool = False,
) -> "pd.DataFrame":
    """
    Compute zonal statistics.

    Parameters
    ----------
    x : SpatRaster
        Values raster.
    z : SpatRaster
        Zones raster (should be integer or categorical).
    fun : str or callable
        Summary function.  Built-ins include ``"sum"``, ``"mean"``, ``"median"``,
        ``"min"``, ``"max"``, ``"sd"``, ``"isNA"`` (count NAs per zone),
        ``"notNA"`` (count non-NA values per zone), …
    na_rm : bool
        Ignore NA values.
    as_raster : bool
        If True, return a SpatRaster instead of a DataFrame.
    wide : bool
        If True (default), return one column per layer of *x*.
    filename : str
    overwrite : bool

    Returns
    -------
    pandas.DataFrame or SpatRaster
    """
    txt = fun if isinstance(fun, str) else getattr(fun, "__name__", "")
    if txt not in _ZONAL_FUNS:
        raise ValueError(f"Function {txt!r} is not supported; use one of {sorted(_ZONAL_FUNS)}")

    opt = spatoptions(filename, overwrite)
    g = _empty_group_raster(x)
    xc = _cpp_zonal(x, z, g, txt, na_rm, opt)
    result = messages(xc, "zonal")

    if as_raster:
        return result

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for zonal()")

    from ._helpers import _getSpatDF
    df = _getSpatDF(result)
    if df is None:
        # result is a SpatRaster; extract via C++ route
        nr, nc_ = result.nrow(), result.ncol()
        import numpy as np
        vals = np.array(result.readValues(0, nr, 0, nc_), dtype=float)
        nms = list(result.names)
        df = pd.DataFrame(vals.reshape(nr * nc_, len(nms)), columns=nms)
    return df
