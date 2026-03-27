"""
zonal.py — zonal statistics.
"""
from __future__ import annotations
from typing import Callable, Optional, Union

from ._terra import SpatRaster, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


_ZONAL_FUNS = {
    "sum", "mean", "median", "modal", "min", "max", "prod", "any", "all",
    "count", "sd", "std", "first",
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
        Summary function.  Built-ins: ``"sum"``, ``"mean"``, ``"median"``,
        ``"min"``, ``"max"``, ``"sd"``, …
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

    opt = SpatOptions(filename, overwrite)
    xc = x.zonal(z, txt, na_rm, opt)
    result = messages(xc, "zonal")

    if as_raster:
        return result

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for zonal()")

    from ._helpers import _getSpatDF
    df = _getSpatDF(result.df) if hasattr(result, 'df') else None
    if df is None:
        # result is a SpatRaster; extract via C++ route
        nr, nc_ = result.nrow(), result.ncol()
        import numpy as np
        vals = np.array(result.readValues(0, nr, 0, nc_), dtype=float)
        nms = list(result.names)
        df = pd.DataFrame(vals.reshape(nr * nc_, len(nms)), columns=nms)
    return df
