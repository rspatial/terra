"""
app.py — apply functions across raster layers.

Covers: app, lapp, tapp, xapp, rapp, sapp.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatOptions
from ._helpers import messages, spatoptions


_CPP_FUNS = {
    "sum", "mean", "median", "modal", "which", "which.min", "which.max",
    "min", "max", "prod", "any", "all", "none", "sd", "std", "first",
}


def _opt() -> SpatOptions:
    return SpatOptions()


def _make_text_fun(fun: Union[str, Callable]) -> Optional[str]:
    """Try to map a callable to its C++ shortcut name."""
    if isinstance(fun, str):
        return fun if fun in _CPP_FUNS else None
    name = getattr(fun, "__name__", "")
    if name in _CPP_FUNS:
        return name
    return None


def app(
    x: SpatRaster,
    fun: Union[str, Callable],
    *,
    na_rm: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a function to each cell across all layers of *x*.

    Built-in string functions (fastest): ``"sum"``, ``"mean"``, ``"median"``,
    ``"min"``, ``"max"``, ``"prod"``, ``"any"``, ``"all"``, ``"none"``,
    ``"sd"``, ``"std"``, ``"first"``, ``"modal"``, ``"which"``,
    ``"which.min"``, ``"which.max"``.

    For any other callable, *fun* receives a 1-D array of length nlyr and
    must return a scalar or 1-D array.

    Parameters
    ----------
    x : SpatRaster
    fun : str or callable
        Aggregation function.
    na_rm : bool
        Ignore NA values.
    filename : str
        Optional output filename.
    overwrite : bool
        Overwrite existing file.

    Returns
    -------
    SpatRaster
    """
    txt = _make_text_fun(fun)
    if txt is not None and txt in _CPP_FUNS:
        opt = spatoptions(filename, overwrite)
        xc = x.summary(txt, na_rm, opt)
        return messages(xc, "app")

    if callable(fun):
        fn = fun
    else:
        raise ValueError(f"Unknown function: {fun!r}")

    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    v = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nl, order='F')
    result = np.apply_along_axis(fn, 1, v)
    if result.ndim == 1:
        result = result.reshape(-1, 1)

    from .rast import rast
    out = rast(x)
    out.nlyr = result.shape[1]
    flat = result.ravel(order='C')
    opt = _opt()
    out.setValues(flat.tolist(), opt)
    return out


def lapp(
    x: SpatRaster,
    fun: Callable,
    *args,
    usenames: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a function to each cell, passing layers as separate arguments.

    The function *fun* is called as ``fun(layer1, layer2, …)`` where each
    argument is a 1-D array of length ncell (one value per cell).

    Parameters
    ----------
    x : SpatRaster
    fun : callable
        Function that accepts nlyr positional array arguments and returns
        an array of length ncell (or ncell × nout).
    *args : additional scalar arguments forwarded to *fun*.
    usenames : bool
        If True, pass layer arrays as keyword arguments named after the
        layer names.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    nms = list(x.names)
    v = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nl, order='F')

    if usenames:
        layer_kwargs = {nms[i]: v[:, i] for i in range(nl)}
        result = fun(*args, **layer_kwargs)
    else:
        layer_args = [v[:, i] for i in range(nl)]
        result = fun(*layer_args, *args)

    result = np.asarray(result, dtype=float)
    if result.ndim == 1:
        result = result.reshape(-1, 1)

    from .rast import rast
    out = rast(x)
    out.nlyr = result.shape[1]
    flat = result.ravel(order='C')
    opt = _opt()
    out.setValues(flat.tolist(), opt)
    return out


def tapp(
    x: SpatRaster,
    index: Union[List, str],
    fun: Union[str, Callable],
    *,
    na_rm: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a function to groups of layers identified by *index*.

    Parameters
    ----------
    x : SpatRaster
    index : list or str
        Group labels of length nlyr, or a time step string such as
        ``"years"``, ``"months"``, ``"days"``, ``"yearmonths"``.
    fun : str or callable
        Aggregation function.
    na_rm : bool
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    nl = x.nlyr()

    # Resolve time-based index strings
    if isinstance(index, str) and x.hasTime:
        import re
        choices = {
            "years": "years", "months": "months", "days": "days",
            "yearmonths": "yearmonths", "yearweeks": "yearweeks",
        }
        step = choices.get(index.lower(), None)
        if step is None:
            raise ValueError(f"Invalid time step: {index!r}")
        raw = list(x.time)
        from datetime import datetime, timezone
        epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
        if step == "years":
            index = [datetime.utcfromtimestamp(t).year for t in raw]
        elif step == "months":
            index = [datetime.utcfromtimestamp(t).month for t in raw]
        elif step == "days":
            index = [datetime.utcfromtimestamp(t).strftime("%Y-%m-%d") for t in raw]
        elif step == "yearmonths":
            index = [datetime.utcfromtimestamp(t).strftime("%Y%m") for t in raw]
        else:
            index = raw

    index = list(index) if not isinstance(index, list) else index
    if len(index) < nl:
        index = (index * ((nl // len(index)) + 1))[:nl]

    txt = _make_text_fun(fun)
    if txt is not None and txt in _CPP_FUNS:
        # Use C++ backend
        int_index = _factor_encode(index)
        nms = _make_unique_names(index)
        opt = spatoptions(filename, overwrite)
        xc = x.apply(int_index, txt, na_rm, nms, [], "", "UTC", opt)
        return messages(xc, "tapp")

    if not callable(fun):
        raise ValueError(f"Unknown function: {fun!r}")

    groups = {}
    for i, g in enumerate(index):
        groups.setdefault(g, []).append(i)

    nr, nc = x.nrow(), x.ncol()
    v_full = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nl, order='F')

    parts = []
    for g in dict.fromkeys(index):
        idx = groups[g]
        chunk = v_full[:, idx]
        r = np.apply_along_axis(fun, 1, chunk)
        if r.ndim == 1:
            r = r.reshape(-1, 1)
        parts.append(r)
    result = np.concatenate(parts, axis=1)

    from .rast import rast
    out = rast(x)
    out.nlyr = result.shape[1]
    opt = _opt()
    out.setValues(result.ravel(order='C').tolist(), opt)
    return out


def xapp(
    x: SpatRaster,
    y: SpatRaster,
    fun: Callable,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a two-argument function cell-by-cell using layers from *x* and *y*.

    *fun* is called as ``fun(x_row, y_row)`` where each argument is a 1-D
    array of length nlyr for one cell.

    Parameters
    ----------
    x : SpatRaster
    y : SpatRaster  Must be compatible with *x*.
    fun : callable
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    nr, nc = x.nrow(), x.ncol()
    nlx = x.nlyr()
    nly = y.nlyr()
    vx = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nlx, order='F')
    vy = np.array(y.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nly, order='F')
    result = np.array([fun(vx[i], vy[i]) for i in range(nr * nc)], dtype=float)
    if result.ndim == 1:
        result = result.reshape(-1, 1)

    from .rast import rast
    out = rast(x)
    out.nlyr = result.shape[1]
    opt = _opt()
    out.setValues(result.ravel(order='C').tolist(), opt)
    return out


def rapp(
    x: SpatRaster,
    first: Union[int, SpatRaster],
    last: Union[int, SpatRaster],
    fun: Union[str, Callable],
    *,
    na_rm: bool = False,
    fill: float = float("nan"),
    clamp: bool = False,
    circular: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a function to a variable-length range of layers per cell.

    Parameters
    ----------
    x : SpatRaster
    first : int or SpatRaster
        First layer index (1-based), or a SpatRaster whose values give
        the start index for each cell.
    last : int or SpatRaster
        Last layer index (1-based), or a SpatRaster.
    fun : str or callable
    na_rm : bool
    fill : float
        Padding value for cells where first > last.
    clamp : bool
        Clamp indices to [1, nlyr].
    circular : bool
        Treat layer sequence as circular.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    if isinstance(first, int) and isinstance(last, int):
        raise ValueError(
            "At least one of 'first' or 'last' must be a SpatRaster. "
            "Use app() for a fixed range."
        )

    firstval = first if isinstance(first, (int, float)) else float("nan")
    lastval = last if isinstance(last, (int, float)) else float("nan")
    index = last if np.isnan(firstval) else first
    if not isinstance(index, SpatRaster):
        raise TypeError("index must be a SpatRaster")

    txt = _make_text_fun(fun)
    if txt and txt in _CPP_FUNS:
        opt = spatoptions(filename, overwrite)
        xc = x.rapply(index, firstval, lastval, txt, clamp, na_rm, circular, opt)
        return messages(xc, "rapp")

    if not callable(fun):
        raise ValueError(f"Unknown function: {fun!r}")

    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    vals = x.rappvals(index, firstval, lastval, clamp, False, fill, 0, nr * nc, circular)
    result = np.array([fun(np.array(v, dtype=float)) for v in vals], dtype=float)
    if result.ndim == 1:
        result = result.reshape(-1, 1)

    from .rast import rast
    out = rast(x)
    out.nlyr = result.shape[1]
    opt = _opt()
    out.setValues(result.ravel(order='C').tolist(), opt)
    return out


def sapp(
    x: SpatRaster,
    fun: Callable,
    *args,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a single-layer function to each layer of *x* separately.

    Parameters
    ----------
    x : SpatRaster
    fun : callable
        Function that accepts a single-layer SpatRaster.
    *args : forwarded to *fun*.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster with the same number of layers as *x*.
    """
    from .rast import rast
    parts = []
    for i in range(x.nlyr()):
        lyr = x.subset([i], _opt())
        parts.append(fun(lyr, *args))
    result = rast(parts)
    return result


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _factor_encode(labels) -> List[int]:
    seen: dict = {}
    result = []
    for v in labels:
        if v not in seen:
            seen[v] = len(seen) + 1
        result.append(seen[v])
    return result


def _make_unique_names(labels) -> List[str]:
    counts: dict = {}
    seen: dict = {}
    result = []
    for v in labels:
        k = str(v)
        if k not in seen:
            seen[k] = 1
            result.append(k)
        else:
            seen[k] += 1
            result.append(f"{k}.{seen[k]}")
    return result
