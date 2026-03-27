"""
math.py — mathematical operations on SpatRaster, SpatExtent, and SpatVector.
"""
from __future__ import annotations
import math as _math
from typing import Optional, Union

from ._terra import SpatRaster, SpatExtent, SpatVector, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# Trigonometric and other single-argument math
# ---------------------------------------------------------------------------

_TRIG_OPS = {
    "sin", "cos", "tan",
    "asin", "acos", "atan",
    "sinh", "cosh", "tanh",
    "asinh", "acosh", "atanh",
    "cospi", "sinpi", "tanpi",
}

_MATH_OPS = {
    "abs", "sign", "sqrt",
    "ceiling", "floor", "trunc",
    "log", "log2", "log10",
    "exp", "expm1", "log1p",
}

_CUM_OPS = {"cumsum", "cumprod", "cummax", "cummin"}


def math(
    x: SpatRaster,
    fun: str,
    digits: int = 0,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Apply a named mathematical function to every cell of *x*.

    Parameters
    ----------
    x : SpatRaster
    fun : str
        One of: ``"abs"``, ``"sign"``, ``"sqrt"``, ``"ceiling"``,
        ``"floor"``, ``"trunc"``, ``"log"``, ``"log2"``, ``"log10"``,
        ``"exp"``, ``"round"``, ``"signif"``,
        ``"sin"``, ``"cos"``, ``"tan"``, ``"asin"``, ``"acos"``,
        ``"atan"``, ``"sinh"``, ``"cosh"``, ``"tanh"``, etc.,
        ``"cumsum"``, ``"cumprod"``, ``"cummax"``, ``"cummin"``.
    digits : int
        Decimal places for ``"round"`` and ``"signif"``.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    opt = SpatOptions(filename, overwrite)
    if fun in _CUM_OPS:
        xc = x.cum(fun[3:], False, opt)
    elif fun in _TRIG_OPS:
        xc = x.trig(fun, opt)
    elif fun in _MATH_OPS:
        xc = x.math(fun, opt)
    elif fun in ("round", "signif"):
        xc = x.math2(fun, int(digits), opt)
    else:
        raise ValueError(f"Unknown math function: {fun!r}")
    return messages(xc, fun)


def log(
    x: SpatRaster,
    base: float = _math.e,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Natural (or custom-base) logarithm.

    Parameters
    ----------
    x : SpatRaster
    base : float
        Log base.  Defaults to e (natural log).
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    opt = SpatOptions(filename, overwrite)
    if base == _math.e:
        xc = x.math("log", opt)
    elif base == 2:
        xc = x.math("log2", opt)
    elif base == 10:
        xc = x.math("log10", opt)
    else:
        import numpy as np
        from .app import app
        return app(x, lambda v: np.log(v) / _math.log(base))
    return messages(xc, "log")


# ---------------------------------------------------------------------------
# Convenience wrappers (one per common operation)
# ---------------------------------------------------------------------------

def sqrt(x: SpatRaster, **kw) -> SpatRaster:
    """Element-wise square root."""
    return math(x, "sqrt", **kw)


def abs_(x: SpatRaster, **kw) -> SpatRaster:
    """Element-wise absolute value."""
    return math(x, "abs", **kw)


def ceiling(x: SpatRaster, **kw) -> SpatRaster:
    """Round up each cell to the nearest integer."""
    return math(x, "ceiling", **kw)


def floor(x: SpatRaster, **kw) -> SpatRaster:
    """Round down each cell to the nearest integer."""
    return math(x, "floor", **kw)


def round_(x: Union[SpatRaster, SpatExtent, SpatVector], digits: int = 0) -> Union[SpatRaster, SpatExtent, SpatVector]:
    """
    Round values to the specified number of decimal places.

    Works for SpatRaster, SpatExtent, and SpatVector.

    Parameters
    ----------
    x : SpatRaster, SpatExtent, or SpatVector
    digits : int

    Returns
    -------
    Same type as *x*.
    """
    if isinstance(x, SpatRaster):
        return math(x, "round", digits=digits)
    elif isinstance(x, SpatExtent):
        xc = x.round(int(digits))
        return xc
    elif isinstance(x, SpatVector):
        xc = x.deepcopy()
        xc.round(int(digits))
        return xc
    raise TypeError(f"round_ does not support {type(x)}")


def cumsum(x: SpatRaster, **kw) -> SpatRaster:
    """Cumulative sum across layers."""
    return math(x, "cumsum", **kw)


def cumprod(x: SpatRaster, **kw) -> SpatRaster:
    """Cumulative product across layers."""
    return math(x, "cumprod", **kw)


def cummax(x: SpatRaster, **kw) -> SpatRaster:
    """Cumulative maximum across layers."""
    return math(x, "cummax", **kw)


def cummin(x: SpatRaster, **kw) -> SpatRaster:
    """Cumulative minimum across layers."""
    return math(x, "cummin", **kw)


# ---------------------------------------------------------------------------
# SpatExtent math
# ---------------------------------------------------------------------------

def floor_ext(x: SpatExtent) -> SpatExtent:
    """Return a copy of *x* with floor-rounded coordinates."""
    return x.floor()


def ceiling_ext(x: SpatExtent) -> SpatExtent:
    """Return a copy of *x* with ceiling-rounded coordinates."""
    return x.ceil()


def round_ext(x: SpatExtent, digits: int = 0) -> SpatExtent:
    """Return a copy of *x* with rounded coordinates."""
    return x.round(int(digits))


# ---------------------------------------------------------------------------
# ifel — raster conditional expression
# ---------------------------------------------------------------------------

def ifel(
    test: SpatRaster,
    yes: Union[SpatRaster, float, int],
    no: Union[SpatRaster, float, int],
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Conditional raster expression: ``where(test, yes, no)``.

    Parameters
    ----------
    test : SpatRaster
        Boolean raster (non-zero = True).
    yes : SpatRaster or scalar
        Values to use where *test* is True.
    no : SpatRaster or scalar
        Values to use where *test* is False.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    from .generics import mask, classify, cover

    yes_is_num = not isinstance(yes, SpatRaster)
    no_is_num = not isinstance(no, SpatRaster)

    if yes_is_num:
        yes_val = float(yes)
    if no_is_num:
        no_val = float(no)

    # Ensure test is 0/1/NA
    logical_test = test.math("abs", _opt())

    if yes_is_num and no_is_num:
        return classify(
            logical_test,
            [[1, yes_val], [0, no_val]],
            filename=filename,
            overwrite=overwrite,
        )

    import math as _m
    if no_is_num:
        no_r = classify(logical_test, [[0, no_val], [1, float("nan")]])
    else:
        no_r = mask(no, logical_test, maskvalues=[float("nan"), 1.0])

    if yes_is_num:
        yes_r = classify(logical_test, [[1, yes_val], [0, float("nan")]])
    else:
        yes_r = mask(yes, logical_test, maskvalues=[float("nan"), 0.0])

    return cover(no_r, yes_r, values=float("nan"), filename=filename, overwrite=overwrite)
