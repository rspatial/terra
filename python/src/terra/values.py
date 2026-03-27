"""
values.py — read/write cell values, metadata, and geometry comparison.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Union, List
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatExtent, SpatOptions

if TYPE_CHECKING:
    pass


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# has_values / in_memory
# ---------------------------------------------------------------------------

def has_values(x: SpatRaster) -> bool:
    """Return True if *x* has cell values loaded."""
    return x.hasValues


def in_memory(x: SpatRaster) -> bool:
    """Return True if *x* is held in memory (not only on disk)."""
    return x.inMemory


def sources(x: SpatRaster) -> List[str]:
    """Return the file path(s) backing *x* (empty string for in-memory)."""
    return list(x.filenames())


# ---------------------------------------------------------------------------
# min / max
# ---------------------------------------------------------------------------

def has_min_max(x: SpatRaster) -> List[bool]:
    """Return per-layer flags indicating whether min/max are stored."""
    return list(x.hasRange)


def min_max(x: SpatRaster, compute: bool = False) -> dict:
    """
    Return the stored or computed min and max values for each layer.

    Parameters
    ----------
    x : SpatRaster
    compute : bool
        If True, compute min/max even if not already stored.

    Returns
    -------
    dict with keys ``"min"`` and ``"max"``, each a list of floats.
    """
    if compute or not all(x.hasRange):
        if compute:
            x.setRange(_opt(), False)
    mins = list(x.range_min)
    maxs = list(x.range_max)
    return {"min": mins, "max": maxs}


def set_min_max(x: SpatRaster, force: bool = False) -> SpatRaster:
    """Compute and store min/max values in *x*."""
    opt = _opt()
    x.setRange(opt, force)
    return x


# ---------------------------------------------------------------------------
# reading values
# ---------------------------------------------------------------------------

def values(
    x: SpatRaster,
    mat: bool = True,
    na_rm: bool = False,
) -> np.ndarray:
    """
    Return all cell values of *x* as a numpy array.

    Parameters
    ----------
    x : SpatRaster
    mat : bool
        If True, return a 2-D matrix (nrows × nlayers); if False a flat 1-D
        vector (all layers concatenated).
    na_rm : bool
        If True, remove rows that contain any NA.

    Returns
    -------
    numpy.ndarray
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    v = x.readValues(0, nr, 0, nc)
    arr = np.array(v, dtype=float)
    if mat:
        arr = arr.reshape(nr * nc, nl, order='F')  # column-major to match R
    if na_rm:
        if arr.ndim == 2:
            mask = ~np.isnan(arr).any(axis=1)
            arr = arr[mask]
        else:
            arr = arr[~np.isnan(arr)]
    return arr


def set_values(x: SpatRaster, v: Union[np.ndarray, List, float]) -> SpatRaster:
    """
    Assign new cell values to a copy of *x*.

    Parameters
    ----------
    x : SpatRaster
        Template raster (geometry is preserved).
    v : array-like or scalar
        Values to assign.  Shape must match ncell × nlyr.

    Returns
    -------
    SpatRaster with new values.
    """
    from .rast import rast
    y = rast(x)
    opt = _opt()
    # pybind11 setValues() only accepts a sequence; broadcast scalars / 0-d arrays.
    if np.ndim(v) == 0:
        val = float(np.asarray(v))
        n = x.nrow() * x.ncol() * x.nlyr()
        y.setValues([val] * n, opt)
    else:
        flat = np.asarray(v, dtype=float).ravel()
        nc = x.nrow() * x.ncol()
        nl = x.nlyr()
        need = nc * nl
        if len(flat) < need:
            flat = np.resize(flat, need)
        elif len(flat) > need:
            flat = flat[:need]
        y.setValues(flat.tolist(), opt)
    return y


# ---------------------------------------------------------------------------
# focal values
# ---------------------------------------------------------------------------

def focal_values(
    x: SpatRaster,
    w: Union[int, List[int]] = 3,
    row: int = 1,
    nrows: int = 1,
    fill: float = float("nan"),
) -> np.ndarray:
    """
    Return focal (neighbourhood) values for a block of rows.

    Parameters
    ----------
    x : SpatRaster
    w : int or list of int
        Window size (cells).  Scalar → square window.
    row : int
        First row (1-based).
    nrows : int
        Number of rows to process.
    fill : float
        Fill value for cells outside the raster.

    Returns
    -------
    numpy.ndarray, shape (nrows × ncol, prod(w)).
    """
    if isinstance(w, (int, float)):
        w = [int(w), int(w)]
    opt = _opt()
    m = x.focalValues(w, fill, max(0, row - 1), nrows, opt)
    arr = np.array(m, dtype=float).reshape(-1, int(w[0]) * int(w[1]))
    return arr


# ---------------------------------------------------------------------------
# SpatVector values
# ---------------------------------------------------------------------------

def vect_values(x: SpatVector) -> "pd.DataFrame":
    """
    Return the attribute table of *x* as a pandas DataFrame.
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for vect_values()")
    from ._helpers import _getSpatDF
    return _getSpatDF(x.df)


def set_vect_values(x: SpatVector, df: "pd.DataFrame") -> SpatVector:
    """
    Set attribute table of a copy of *x*.

    Parameters
    ----------
    x : SpatVector
    df : pandas.DataFrame
        New attribute table.  Number of rows must match nrow(x).
    """
    from ._helpers import _makeSpatDF
    xc = x.deepcopy()
    sdf = _makeSpatDF(df)
    xc.set_df(sdf)
    return xc


# ---------------------------------------------------------------------------
# compare_geom
# ---------------------------------------------------------------------------

def compare_geom(
    x: SpatRaster,
    y: SpatRaster,
    lyrs: bool = False,
    crs: bool = True,
    warncrs: bool = False,
    ext: bool = True,
    rowcol: bool = True,
    res: bool = False,
    stop_on_error: bool = True,
    tolerance: Optional[float] = None,
) -> bool:
    """
    Compare the geometry of two SpatRasters.

    Parameters
    ----------
    x, y : SpatRaster
    lyrs : bool
        Check number of layers.
    crs : bool
        Check CRS.
    warncrs : bool
        Warn instead of error on CRS mismatch.
    ext : bool
        Check extent.
    rowcol : bool
        Check number of rows and columns.
    res : bool
        Check resolution.
    stop_on_error : bool
        Raise RuntimeError on failure.
    tolerance : float, optional
        Coordinate tolerance.

    Returns
    -------
    bool
    """
    if tolerance is None:
        opt = spatoptions()
        tolerance = opt.tolerance
    result = x.compare_geom(y, lyrs, crs, tolerance, warncrs, ext, rowcol, res)
    if not result and stop_on_error:
        msgs = []
        if x.has_error():
            msgs.append(x.getError())
        raise RuntimeError("; ".join(msgs) or "Geometries do not match")
    return bool(result)
