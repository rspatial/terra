"""
merge.py — merge and mosaic rasters and join vector attribute tables.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union

from ._terra import SpatRaster, SpatRasterCollection, SpatVector, SpatOptions
from ._helpers import messages, spatoptions


def _opt() -> SpatOptions:
    return SpatOptions()


def _sprc_from_rasters(rasters: List[SpatRaster]) -> SpatRasterCollection:
    """Build a :class:`SpatRasterCollection` from rasters (R ``sprc(...)``)."""
    rc = SpatRasterCollection()
    for r in rasters:
        rc.add(r.deepcopy(), "")
    return rc


# ---------------------------------------------------------------------------
# Raster merge / mosaic
# ---------------------------------------------------------------------------

def merge(
    x: SpatRaster,
    *others: SpatRaster,
    first: bool = True,
    na_rm: bool = True,
    algo: int = 1,
    method: Optional[str] = None,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Merge two or more SpatRasters that overlap or are adjacent.

    Where rasters overlap the value from *x* is used when *first=True*
    (default) or the last raster with a non-NA value when *first=False*.

    Parameters
    ----------
    x : SpatRaster
    *others : SpatRaster
        Additional rasters to merge.
    first : bool
        Use the first raster's value in overlapping areas.
    na_rm : bool
        Skip NA values when merging.
    algo : int
        Internal algorithm (1 or 2).
    method : str, optional
        Blending method for smooth edges.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    all_rasters = [x] + list(others)
    opt = spatoptions(filename, overwrite)

    rc = _sprc_from_rasters(all_rasters)
    if method is None:
        method = ""
    xc = rc.merge(first, na_rm, algo, method, opt)
    return messages(xc, "merge")


def mosaic(
    x: SpatRaster,
    *others: SpatRaster,
    fun: Union[str, Callable] = "mean",
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Mosaic two or more overlapping SpatRasters using an aggregation function.

    Unlike merge(), mosaic() aggregates overlapping values rather than
    giving priority to any one raster.

    Parameters
    ----------
    x : SpatRaster
    *others : SpatRaster
    fun : str or callable
        Aggregation function: ``"mean"`` (default), ``"sum"``, ``"min"``,
        ``"max"``, ``"first"``, ``"last"``.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    fun_str = fun if isinstance(fun, str) else getattr(fun, "__name__", "mean")
    if fun_str not in {"mean", "sum", "min", "max", "first", "last"}:
        raise ValueError(f"fun must be one of mean/sum/min/max/first/last; got {fun_str!r}")

    all_rasters = [x] + list(others)
    opt = spatoptions(filename, overwrite)
    rc = _sprc_from_rasters(all_rasters)
    xc = rc.mosaic(fun_str, opt)
    return messages(xc, "mosaic")


# ---------------------------------------------------------------------------
# Vector attribute table join
# ---------------------------------------------------------------------------

def merge_vect(
    x: SpatVector,
    y: "pd.DataFrame",
    **kwargs,
) -> SpatVector:
    """
    Join attribute columns from *y* to *x*.

    Parameters
    ----------
    x : SpatVector
    y : pandas.DataFrame
        Table to join.  Must share at least one column with *x*'s
        attribute table.
    **kwargs
        Forwarded to ``pandas.merge`` (e.g. ``on``, ``how``).

    Returns
    -------
    SpatVector
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for merge_vect()")
    from ._helpers import _getSpatDF, _makeSpatDF

    v = _getSpatDF(x.df)
    if v is None:
        v = pd.DataFrame()
    uid_col = "__uid__"
    v[uid_col] = range(len(v))
    m = pd.merge(v, y, **kwargs)
    m = m.sort_values(uid_col).reset_index(drop=True)
    if len(m) > len(x):
        raise ValueError("merge would expand the number of features; 'all.y=True' is not supported")
    uid_vals = m[uid_col].dropna().astype(int).tolist()
    m = m.drop(columns=[uid_col])
    xc = x.subset(uid_vals, _opt())
    sdf = _makeSpatDF(m)
    xc.set_df(sdf)
    return xc
