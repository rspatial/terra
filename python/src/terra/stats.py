"""
stats.py — row/column statistics and value matching for SpatRaster.

Covers functionality from rowSums.R, match.R, autocor.R, layerCor.R.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatOptions
from ._helpers import messages, spatoptions


def _opt() -> SpatOptions:
    return SpatOptions()


def _read_values_layer_matrix(x: SpatRaster) -> np.ndarray:
    """
    Read all cell values as (nrow * ncol, nlyr), Fortran order (aligned with R).

    GDAL-backed rasters require readStart before readValues.
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    x.readStart()
    try:
        raw = x.readValues(0, nr, 0, nc)
    finally:
        x.readStop()
    return np.array(raw, dtype=float).reshape(nr * nc, nl, order="F")


# ---------------------------------------------------------------------------
# Row / column statistics
# ---------------------------------------------------------------------------

def row_sums(x: SpatRaster, na_rm: bool = False) -> np.ndarray:
    """
    Sum of each row across columns, returned separately per layer.

    Parameters
    ----------
    x : SpatRaster
    na_rm : bool

    Returns
    -------
    numpy.ndarray, shape (nrow, nlyr).
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    vals = _read_values_layer_matrix(x)
    vals_3d = vals.reshape(nr, nc, nl)
    if na_rm:
        return np.nansum(vals_3d, axis=1)
    return np.sum(vals_3d, axis=1)


def col_sums(x: SpatRaster, na_rm: bool = False) -> np.ndarray:
    """
    Sum of each column across rows, returned separately per layer.

    Parameters
    ----------
    x : SpatRaster
    na_rm : bool

    Returns
    -------
    numpy.ndarray, shape (ncol, nlyr).
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    vals = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nl, order='F')
    vals_3d = vals.reshape(nr, nc, nl)
    if na_rm:
        return np.nansum(vals_3d, axis=0)
    return np.sum(vals_3d, axis=0)


def row_means(x: SpatRaster, na_rm: bool = False) -> np.ndarray:
    """
    Mean of each row across columns, returned separately per layer.

    Parameters
    ----------
    x : SpatRaster
    na_rm : bool

    Returns
    -------
    numpy.ndarray, shape (nrow, nlyr).
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    vals = _read_values_layer_matrix(x)
    vals_3d = vals.reshape(nr, nc, nl)
    if na_rm:
        return np.nanmean(vals_3d, axis=1)
    return np.mean(vals_3d, axis=1)


def col_means(x: SpatRaster, na_rm: bool = False) -> np.ndarray:
    """
    Mean of each column across rows, returned separately per layer.

    Parameters
    ----------
    x : SpatRaster
    na_rm : bool

    Returns
    -------
    numpy.ndarray, shape (ncol, nlyr).
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    vals = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nl, order='F')
    vals_3d = vals.reshape(nr, nc, nl)
    if na_rm:
        return np.nanmean(vals_3d, axis=0)
    return np.mean(vals_3d, axis=0)


# ---------------------------------------------------------------------------
# match / is_in
# ---------------------------------------------------------------------------

def match_rast(
    x: SpatRaster,
    table: List,
    nomatch: float = float("nan"),
) -> SpatRaster:
    """
    Return the position of each cell's value in *table*.

    Parameters
    ----------
    x : SpatRaster
    table : list
        Values to match against.
    nomatch : float
        Value to use when there is no match (default: NaN).

    Returns
    -------
    SpatRaster
    """
    import numpy as np
    from .app import app

    table_u = list(dict.fromkeys(table))
    table_np = np.array(table_u, dtype=float)

    def _match(v: np.ndarray) -> np.ndarray:
        out = np.full_like(v, nomatch, dtype=float)
        for idx, t in enumerate(table_np):
            mask = v == t
            out[mask] = idx + 1
        return out

    return app(x, _match)


def is_in(
    x: SpatRaster,
    table: List,
) -> SpatRaster:
    """
    Return a binary raster: 1 where cell values are in *table*, 0 otherwise.

    Parameters
    ----------
    x : SpatRaster
    table : list of values

    Returns
    -------
    SpatRaster
    """
    table_u = list(dict.fromkeys([t for t in table if t is not None]))
    opt = _opt()
    xc = x.is_in([float(t) for t in table_u], opt)
    return messages(xc, "is_in")


# ---------------------------------------------------------------------------
# autocor — spatial autocorrelation (Moran's I)
# ---------------------------------------------------------------------------

def autocor(
    x: SpatRaster,
    w: Optional[Union[str, "np.ndarray"]] = None,
    global_: bool = True,
    method: str = "moran",
    filename: str = "",
    overwrite: bool = False,
) -> Union[float, SpatRaster]:
    """
    Compute spatial autocorrelation for *x*.

    Parameters
    ----------
    x : SpatRaster
        Single-layer raster.
    w : str or array, optional
        Weights matrix.  ``"queen"`` (8-neighbour, default), ``"rook"``
        (4-neighbour), or a custom 2-D weight matrix.
    global_ : bool
        If True, return a scalar (global Moran's I or Geary's C).
        If False, return a local SpatRaster.
    method : str
        ``"moran"`` (default) or ``"geary"``.
    filename : str
    overwrite : bool

    Returns
    -------
    float (global) or SpatRaster (local).
    """
    if w is None:
        w = "queen"
    opt = spatoptions(filename, overwrite)
    xc = x.autocor(w if isinstance(w, str) else "custom", not global_, method, opt)
    result = messages(xc, "autocor")
    if global_:
        # Returns a single-cell raster; extract the value
        result.readStart()
        try:
            v = result.readValues(0, 1, 0, 1)
        finally:
            result.readStop()
        return float(v[0]) if v else float("nan")
    return result


# ---------------------------------------------------------------------------
# layerCor — layer-wise correlation matrix
# ---------------------------------------------------------------------------

def layer_cor(
    x: SpatRaster,
    fun: str = "pearson",
    *,
    na_rm: bool = True,
    asSample: bool = True,
) -> "np.ndarray":
    """
    Compute the correlation (or covariance) matrix between layers of *x*.

    Parameters
    ----------
    x : SpatRaster
    fun : str
        ``"pearson"`` (default), ``"spearman"``, or ``"cov"``
        (covariance).
    na_rm : bool
        Ignore NA values.
    asSample : bool
        Use sample (n-1) divisor for covariance.

    Returns
    -------
    numpy.ndarray, shape (nlyr, nlyr).
    """
    nl = x.nlyr()
    vals = _read_values_layer_matrix(x)
    if na_rm:
        mask = ~np.isnan(vals).any(axis=1)
        vals = vals[mask]

    if fun == "cov":
        ddof = 1 if asSample else 0
        return np.cov(vals.T, ddof=ddof)
    elif fun == "spearman":
        from scipy.stats import spearmanr
        corr, _ = spearmanr(vals)
        if nl == 1:
            return np.array([[1.0]])
        return np.array(corr)
    else:
        return np.corrcoef(vals.T)
