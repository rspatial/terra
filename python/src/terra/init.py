"""
init.py — initialise SpatRaster cells with coordinate or positional values.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


_INIT_STRINGS = {"x", "y", "xy", "row", "col", "cell", "chess"}


def init(
    x: SpatRaster,
    fun: Union[str, float, int, bool, np.ndarray, List, Callable],
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Initialise a SpatRaster with a coordinate-based or constant value.

    Parameters
    ----------
    x : SpatRaster
        Template raster whose geometry is used.
    fun : str, scalar, array-like, or callable
        - ``"x"`` — cell x-coordinate.
        - ``"y"`` — cell y-coordinate.
        - ``"xy"`` — both (returns 2-layer raster).
        - ``"row"`` — row number (1-based).
        - ``"col"`` — column number (1-based).
        - ``"cell"`` — cell number (1-based).
        - ``"chess"`` — checkerboard pattern (0/1).
        - scalar — fill all cells with that value.
        - array-like — fill cells with given values (recycled as needed).
        - callable — called as ``fun(n)`` where n = ncell * nlyr, must
          return an array of that length.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    from .rast import rast

    opt = SpatOptions(filename, overwrite)

    if isinstance(fun, str):
        fun = fun[0] if len(fun) == 0 else fun
        if fun not in _INIT_STRINGS:
            raise ValueError(
                f"Unknown init function {fun!r}. "
                f"Choose from {sorted(_INIT_STRINGS)}."
            )
        y = rast(x, nlyrs=1)
        xc = y.initf(fun, True, opt)
        return messages(xc, "init")

    if callable(fun):
        nc = x.ncol() * x.nlyr()
        b_size = x.nrow() * nc
        result = fun(b_size)
        result = np.asarray(result, dtype=float)
        y = rast(x)
        y.setValues(result.tolist(), opt)
        return y

    # Scalar or array-like
    y = rast(x)
    if np.isscalar(fun):
        val = float(fun)
        xc = y.initv(val, opt)
        return messages(xc, "init")

    arr = np.asarray(fun, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == x.ncol():
        arr = arr.ravel(order='C')
    total = x.nrow() * x.ncol() * x.nlyr()
    if len(arr) < total:
        arr = np.resize(arr, total)
    elif len(arr) > total:
        arr = arr[:total]
    xc = y.initv(arr.tolist(), opt)
    return messages(xc, "init")
