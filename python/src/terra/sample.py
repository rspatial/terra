"""
sample.py — spatial sampling of SpatRaster and SpatExtent.
"""
from __future__ import annotations
from typing import Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatExtent, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


def spat_sample(
    x: Union[SpatRaster, SpatExtent],
    size: int,
    method: str = "random",
    *,
    replace: bool = False,
    na_rm: bool = False,
    as_raster: bool = False,
    as_points: bool = False,
    values: bool = True,
    cells: bool = False,
    xy: bool = False,
    lonlat: Optional[bool] = None,
    exact: bool = False,
    warn: bool = True,
) -> Union["pd.DataFrame", SpatVector, SpatRaster, np.ndarray]:
    """
    Draw a spatial sample from a SpatRaster or SpatExtent.

    Parameters
    ----------
    x : SpatRaster or SpatExtent
    size : int
        Number of samples.
    method : str
        ``"random"`` (default), ``"regular"``, or ``"stratified"``.
    replace : bool
        Sample with replacement.
    na_rm : bool
        Exclude NA cells.
    as_raster : bool
        Return sampled cells as a SpatRaster mask.
    as_points : bool
        Return sampled cells as a SpatVector of points.
    values : bool
        Include cell values in the output.
    cells : bool
        Include cell numbers.
    xy : bool
        Include x/y coordinates.
    lonlat : bool, optional
        Override CRS-based detection for geographic correction.
    exact : bool
        Use exact regular grid (SpatExtent only).
    warn : bool
        Warn when fewer samples than requested are available.

    Returns
    -------
    pandas.DataFrame, SpatVector, SpatRaster, or numpy.ndarray
    """
    size = max(1, int(round(size)))

    if isinstance(x, SpatExtent):
        return _sample_extent(x, size, method, lonlat, as_points, exact)

    method = method.lower()
    if method not in ("random", "regular", "stratified"):
        raise ValueError(f"method must be 'random', 'regular', or 'stratified'; got {method!r}")

    if lonlat is None:
        lonlat = x.is_lonlat() if hasattr(x, "is_lonlat") else False

    opt = _opt()
    result = x.spatSample(size, method, replace, na_rm, as_raster, as_points,
                          False, xy, cells, lonlat, exact, warn, opt)

    if as_raster or as_points:
        return messages(result, "spatSample")

    # Build DataFrame
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for spat_sample()")
    from ._helpers import _getSpatDF
    df = _getSpatDF(result)
    if df is None:
        return pd.DataFrame()
    return df


def _sample_extent(
    x: SpatExtent,
    size: int,
    method: str,
    lonlat: Optional[bool],
    as_points: bool,
    exact: bool,
) -> Union[np.ndarray, SpatVector]:
    method = method.lower()
    if method not in ("random", "regular"):
        raise ValueError("method for SpatExtent must be 'random' or 'regular'")
    if lonlat is None:
        raise ValueError("lonlat must be specified when sampling from a SpatExtent")
    if method == "random":
        s = x.sampleRandom(size, lonlat, 0)
    else:
        s = x.sampleRegular(size, lonlat)
    arr = np.array(s, dtype=float).reshape(-1, 2)
    if as_points:
        from .vect import vect
        return vect(arr, type="points")
    return arr


# ---------------------------------------------------------------------------
# grid_sample — spatial thinning on a grid
# ---------------------------------------------------------------------------

def grid_sample(
    xy: Union[np.ndarray, SpatVector],
    r: SpatRaster,
    n: int = 1,
    chess: str = "",
) -> np.ndarray:
    """
    Thin point locations to at most *n* per raster cell.

    Parameters
    ----------
    xy : ndarray (n, 2) or SpatVector of points
        Input coordinates.
    r : SpatRaster
        Reference raster defining the grid.
    n : int
        Maximum number of points per cell.
    chess : str
        ``""`` (all cells), ``"white"`` or ``"black"`` (checkerboard
        sub-selection).

    Returns
    -------
    numpy.ndarray of indices (1-based) into *xy*.
    """
    if isinstance(xy, SpatVector):
        crds = np.array(xy.coordinates(), dtype=float).reshape(-1, 2)
    else:
        crds = np.asarray(xy, dtype=float)
        if crds.ndim == 1:
            crds = crds.reshape(1, -1)

    cell = r.cellFromXY(crds[:, 0].tolist(), crds[:, 1].tolist(), float("nan"))
    cell = np.array(cell, dtype=float) + 1
    valid = ~np.isnan(cell)
    cell = cell.astype(int)

    # Optional checkerboard filter
    if chess.lower() in ("white", "black"):
        nc = r.ncol()
        row_idx = (cell - 1) // nc
        col_idx = (cell - 1) % nc
        parity = (row_idx + col_idx) % 2
        keep_parity = 0 if chess.lower() == "white" else 1
        valid &= (parity == keep_parity)

    selected = []
    from collections import defaultdict
    cell_pts: dict = defaultdict(list)
    for i, (c, v) in enumerate(zip(cell, valid)):
        if v:
            cell_pts[c].append(i)

    rng = np.random.default_rng()
    for c, idxs in cell_pts.items():
        if len(idxs) <= n:
            selected.extend(idxs)
        else:
            chosen = rng.choice(idxs, n, replace=False)
            selected.extend(chosen.tolist())

    return np.array(sorted(selected), dtype=int) + 1
