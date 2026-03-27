"""
focal.py — focal (moving-window) statistics and focal matrices.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatOptions
from ._helpers import messages

_cpp_focal = SpatRaster.focal  # captured before monkey-patching


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# focal matrix helpers
# ---------------------------------------------------------------------------

def focal_mat(
    type_: str = "circle",
    r: Union[float, List[float]] = 1,
    lp: float = 2.0,
) -> np.ndarray:
    """
    Create a focal weight matrix.

    Parameters
    ----------
    type_ : str
        One of ``"circle"``, ``"square"``, ``"rectangle"``, ``"Gauss"``,
        ``"band"``.
    r : float or list of float
        Radius or size parameters (interpretation depends on *type_*).
    lp : float
        Exponent for ``"band"`` type.

    Returns
    -------
    numpy.ndarray  (dtype float64, NaN = outside window)
    """
    if isinstance(r, (int, float)):
        r = [float(r)]
    else:
        r = [float(v) for v in r]

    t = type_.lower()
    if t in ("circle", "disc"):
        rad = r[0]
        size = int(2 * rad + 1)
        arr = np.full((size, size), float("nan"))
        cx, cy = rad, rad
        for i in range(size):
            for j in range(size):
                if (i - cx) ** 2 + (j - cy) ** 2 <= rad ** 2:
                    arr[i, j] = 1.0
        return arr
    elif t == "square":
        s = int(r[0])
        if s % 2 == 0:
            s += 1
        return np.ones((s, s), dtype=float)
    elif t == "rectangle":
        if len(r) < 2:
            r = [r[0], r[0]]
        nr_, nc_ = int(r[0]), int(r[1])
        if nr_ % 2 == 0:
            nr_ += 1
        if nc_ % 2 == 0:
            nc_ += 1
        return np.ones((nr_, nc_), dtype=float)
    elif t in ("gauss", "gaussian"):
        sigma = r[0]
        s = int(2 * np.ceil(3 * sigma) + 1)
        if s % 2 == 0:
            s += 1
        c = s // 2
        arr = np.zeros((s, s), dtype=float)
        for i in range(s):
            for j in range(s):
                arr[i, j] = np.exp(-((i - c) ** 2 + (j - c) ** 2) / (2 * sigma ** 2))
        return arr / arr.sum()
    elif t == "band":
        min_d = r[0] if len(r) > 0 else 0
        max_d = r[1] if len(r) > 1 else r[0]
        s = int(2 * max_d + 1)
        arr = np.full((s, s), float("nan"))
        c = s // 2
        for i in range(s):
            for j in range(s):
                d = np.sqrt((i - c) ** 2 + (j - c) ** 2)
                if min_d <= d <= max_d:
                    arr[i, j] = 1.0
        return arr
    else:
        raise ValueError(f"Unknown focal matrix type: {type_!r}")


# ---------------------------------------------------------------------------
# focal
# ---------------------------------------------------------------------------

_FOCAL_FUNS = {
    "sum", "mean", "median", "modal", "min", "max", "prod", "any", "all",
    "count", "sd", "std", "first",
}


def focal(
    x: SpatRaster,
    w: Union[int, List, np.ndarray],
    fun: Union[str, Callable] = "sum",
    *,
    na_rm: bool = True,
    fillvalue: float = float("nan"),
    expand: bool = False,
    na_policy: str = "omit",
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a focal (moving-window) function to *x*.

    Parameters
    ----------
    x : SpatRaster
    w : int, list, or array
        Window definition.  An integer is treated as a square window of
        that size; a 2-element list [rows, cols] defines a rectangular
        window; a 2-D array defines a weighted focal filter.
    fun : str or callable
        Aggregation function.  Built-in shortcuts: ``"sum"``, ``"mean"``,
        ``"median"``, ``"min"``, ``"max"``, ``"sd"``, ``"modal"``, …
    na_rm : bool
        Ignore NA values.
    fillvalue : float
        Value to use for cells outside the raster border.
    expand : bool
        Expand the output to accommodate window overlap.
    na_policy : str
        ``"omit"`` (skip NA input), ``"only"`` (only apply to NA cells),
        or ``"all"`` (apply even when any NA is present).
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    if isinstance(w, (int, float)):
        wmat = [int(w), int(w)]
        wvals = [1.0] * (int(w) ** 2)
    elif isinstance(w, (list, tuple)) and all(isinstance(v, (int, float)) for v in w) and len(w) <= 2:
        r, c = (int(w[0]), int(w[1])) if len(w) == 2 else (int(w[0]), int(w[0]))
        wmat = [r, c]
        wvals = [1.0] * (r * c)
    else:
        arr = np.asarray(w, dtype=float)
        if arr.ndim != 2:
            raise ValueError("w must be a scalar, 2-element list, or 2-D array")
        wmat = list(arr.shape)
        wvals = arr.ravel(order='C').tolist()

    txt = fun if isinstance(fun, str) else getattr(fun, "__name__", "")
    if txt in _FOCAL_FUNS:
        opt = SpatOptions(filename, overwrite)
        xc = _cpp_focal(x, wmat, wvals, fillvalue, na_rm, txt, expand, na_policy, opt)
        return messages(xc, "focal")

    if not callable(fun):
        raise ValueError(f"Unknown focal function: {fun!r}")

    opt = _opt()
    fv = x.focalValues(wmat, fillvalue, 0, x.nrow(), opt)
    nr, nc_ = x.nrow(), x.ncol()
    n_w = wmat[0] * wmat[1]
    fv_np = np.array(fv, dtype=float).reshape(nr * nc_, n_w)
    result = np.array([fun(fv_np[i]) for i in range(nr * nc_)], dtype=float)

    from .rast import rast
    out = rast(x)
    out.setValues(result.tolist(), opt)
    return out


# ---------------------------------------------------------------------------
# focal3D
# ---------------------------------------------------------------------------

def focal3D(
    x: SpatRaster,
    w: List,
    fun: Union[str, Callable] = "sum",
    *,
    na_rm: bool = True,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    3-D focal function across rows, columns, and layers.

    Parameters
    ----------
    x : SpatRaster
    w : list [nlyr_window, nrow_window, ncol_window]
    fun : str or callable
    na_rm : bool
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    if len(w) != 3:
        raise ValueError("w must have three elements [nlyr, nrow, ncol]")
    txt = fun if isinstance(fun, str) else getattr(fun, "__name__", "")
    if txt in _FOCAL_FUNS:
        opt = SpatOptions(filename, overwrite)
        xc = x.focal3d(w, txt, na_rm, opt)
        return messages(xc, "focal3D")
    raise NotImplementedError("Custom functions are not supported for focal3D")
