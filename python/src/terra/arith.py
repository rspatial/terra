"""
Arithmetic, comparison, logic, and NA operators for terra objects.

Mirrors R ``Arith_generics.R``.

**SpatExtent operators** (registered via :func:`register_operators`):
  ``e + n``  expand   ``e - n``  shrink   ``e * n``  scale from centre
  ``e / n``  shrink   ``e % n``  align    ``e + e2`` union
  ``e * e2`` intersect  ``e / e2``  returns (dx_ratio, dy_ratio)
  ``e == e2`` etc.  comparison

**SpatRaster operators**:
  ``r + r2 / n``  ``r - …``  ``r * …``  ``r / …``  ``r ** …``
  ``r == …``  ``r != …``  ``r > …``  ``r < …``  ``r >= …``  ``r <= …``
  ``r & …``  ``r | …``  ``~r``  ``-r``

**SpatVector operators**:
  ``v + v2``  union   ``v * v2``  intersect   ``v - v2``  erase

Standalone functions (also exported from ``terra``):
  :func:`is_na`, :func:`not_na`, :func:`is_true`, :func:`is_false`,
  :func:`is_nan`, :func:`is_finite`, :func:`is_infinite`,
  :func:`any_na`, :func:`all_na`, :func:`no_na`, :func:`count_na`,
  :func:`which_max`, :func:`which_min`, :func:`which_lyr`,
  :func:`where_max`, :func:`where_min`,
  :func:`rast_sum`, :func:`rast_mean`, :func:`rast_min`, :func:`rast_max`,
  :func:`rast_median`, :func:`rast_modal`,
  :func:`compare_rast`, :func:`logic_rast`,
  :func:`as_int_rast`, :func:`as_bool_rast`,
  :func:`is_bool_rast`, :func:`is_int_rast`, :func:`is_num_rast`,
  :func:`stdev_rast`,
"""

from __future__ import annotations

import math
from typing import Any, List, Optional, Union

from ._helpers import messages
from ._terra import SpatExtent, SpatOptions, SpatRaster, SpatVector

__all__ = [
    # NA / logical tests
    "is_na", "not_na", "is_true", "is_false",
    "is_nan", "is_finite", "is_infinite",
    "any_na", "all_na", "no_na", "count_na",
    # summaries
    "which_max", "which_min", "which_lyr",
    "where_max", "where_min",
    "rast_sum", "rast_mean", "rast_min", "rast_max",
    "rast_median", "rast_modal",
    "stdev_rast",
    # compare / logic
    "compare_rast", "logic_rast_fn",
    # type coercion / queries
    "as_int_rast", "as_bool_rast",
    "is_bool_rast", "is_int_rast", "is_num_rast",
    # operator registration
    "register_operators",
]


# ── Helpers ──────────────────────────────────────────────────────────────────

def _opt() -> SpatOptions:
    return SpatOptions()


def _to_floats(x: Any) -> Union[float, List[float]]:
    """Accept scalar or sequence, return float or list[float]."""
    if isinstance(x, (int, float)):
        return [float(x)]
    return [float(v) for v in x]


# ── SpatExtent arithmetic ────────────────────────────────────────────────────

def _ext_add(e: SpatExtent, n: Any) -> SpatExtent:
    """e + n  →  expand (grow each side by n)."""
    from .extent import ext as _ext
    if isinstance(n, SpatExtent):
        # union
        ec = e.deepcopy()
        ec.union(n)
        return ec
    v = e.vector
    ns = [float(x) for x in _to_floats(n)]
    while len(ns) < 4:
        ns.extend(ns)
    ns = ns[:4]
    # R: e2[c(1,3)] <- -e2[c(1,3)]; then e + e2
    # meaning xmin -= n[0], xmax += n[1], ymin -= n[2], ymax += n[3]
    return _ext(v[0] - ns[0], v[1] + ns[1], v[2] - ns[2], v[3] + ns[3])


def _ext_sub(e: SpatExtent, n: Any) -> SpatExtent:
    """e - n  →  shrink."""
    from .extent import ext as _ext
    v = e.vector
    ns = [float(x) for x in _to_floats(n)]
    while len(ns) < 4:
        ns.extend(ns)
    ns = ns[:4]
    result = _ext(v[0] + ns[0], v[1] - ns[1], v[2] + ns[2], v[3] - ns[3])
    if not result.valid:
        raise ValueError("- : this would create an invalid extent")
    return result


def _ext_mul(e: SpatExtent, n: Any) -> Any:
    """e * n  →  scale from centre; e * e2 → intersect."""
    from .extent import ext as _ext
    if isinstance(n, SpatExtent):
        result = e.intersect(n)
        if not result.valid_notempty:
            return None
        return result
    ns = [abs(float(x)) for x in _to_floats(n)]
    while len(ns) < 4:
        ns.extend(ns)
    ns = ns[:4]
    v = e.vector
    dx = v[1] - v[0]
    dy = v[3] - v[2]
    mx = v[0] + dx / 2
    my = v[2] + dy / 2
    result = _ext(mx - (dx / 2) * ns[0], mx + (dx / 2) * ns[1],
                  my - (dy / 2) * ns[2], my + (dy / 2) * ns[3])
    if not result.valid:
        raise ValueError("* : this would create an invalid extent")
    return result


def _ext_div(e: SpatExtent, n: Any) -> Any:
    """e / n  →  shrink from centre; e / e2 → (dx_ratio, dy_ratio)."""
    from .extent import ext as _ext
    if isinstance(n, SpatExtent):
        v1, v2 = e.vector, n.vector
        return ((v1[1] - v1[0]) / (v2[1] - v2[0]),
                (v1[3] - v1[2]) / (v2[3] - v2[2]))
    ns = [abs(float(x)) for x in _to_floats(n)]
    while len(ns) < 4:
        ns.extend(ns)
    ns = ns[:4]
    v = e.vector
    dx = v[1] - v[0]
    dy = v[3] - v[2]
    mx = v[0] + dx / 2
    my = v[2] + dy / 2
    result = _ext(mx - dx / (2 * ns[0]), mx + dx / (2 * ns[1]),
                  my - dy / (2 * ns[2]), my + dy / (2 * ns[3]))
    if not result.valid:
        raise ValueError("/ : this would create an invalid extent")
    return result


def _ext_mod(e: SpatExtent, n: Any) -> SpatExtent:
    """e % n  →  align to resolution n."""
    return e.align(float(n), "")


def _ext_compare(e: SpatExtent, other: SpatExtent, op: str) -> bool:
    if not isinstance(other, SpatExtent):
        return NotImplemented
    return e.compare(other, op, 1e-6)


# ── SpatRaster arithmetic ─────────────────────────────────────────────────────

def _rast_arith_rast(a: SpatRaster, b: SpatRaster, op: str) -> SpatRaster:
    r = a.arith_rast(b, op, False, _opt())
    return messages(r, op)


def _rast_arith_numb(a: SpatRaster, b: Any, op: str, reverse: bool = False) -> SpatRaster:
    vals = _to_floats(b)
    r = a.arith_numb(vals, op, reverse, False, _opt())
    return messages(r, op)


def _rast_binop(a: SpatRaster, b: Any, op: str) -> SpatRaster:
    if isinstance(b, SpatRaster):
        return _rast_arith_rast(a, b, op)
    try:
        import numpy as np  # type: ignore
        if isinstance(b, np.ndarray):
            flat = b.astype("float64").ravel(order="C").tolist()
            dims = [int(b.shape[0]), int(b.shape[1])] if b.ndim == 2 else [1, int(b.size)]
            r = a.arith_m(flat, op, dims, False, _opt())
            return messages(r, op)
    except ImportError:
        pass
    if isinstance(b, bool):
        b = int(b)
    return _rast_arith_numb(a, b, op, False)


def _rast_rbinop(b: Any, a: SpatRaster, op: str) -> SpatRaster:
    """numeric OP raster  (reversed)."""
    if isinstance(b, bool):
        b = int(b)
    return _rast_arith_numb(a, b, op, True)


def _rast_logic_rast(a: SpatRaster, b: SpatRaster, op: str) -> SpatRaster:
    r = a.logic_rast(b, op, _opt())
    return messages(r, op)


def _rast_logic_numb(a: SpatRaster, b: Any, op: str) -> SpatRaster:
    r = a.logic_numb([bool(b)], op, _opt())
    return messages(r, op)


def _rast_logic_op(a: SpatRaster, b: Any, op: str) -> SpatRaster:
    if isinstance(b, SpatRaster):
        return _rast_logic_rast(a, b, op)
    return _rast_logic_numb(a, bool(b), op)


# ── SpatVector arithmetic ────────────────────────────────────────────────────

def _vect_binop(a: SpatVector, b: SpatVector, op: str) -> SpatVector:
    if a.type() != b.type():
        raise TypeError(f"{op}: geometry types do not match")
    if op == "+":
        r = a.union(b)
    elif op == "*":
        r = a.intersect(b, True)
    elif op == "-":
        r = a.erase_agg(b)
    else:
        raise TypeError(f"{op}: only +, *, - are supported for SpatVector")
    return messages(r, op)


# ── NA / logical tests ────────────────────────────────────────────────────────

def _wrap(r: SpatRaster, name: str, false_na: bool = False) -> SpatRaster:
    r = r.is_wrapper(name, false_na, _opt())
    return messages(r, name)


def is_na(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Cell-wise ``is.na()`` — like R ``is.na(SpatRaster)``."""
    return _wrap(x, "isnan")


def not_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """Cell-wise ``not.na()``."""
    return _wrap(x, "isnotnan", false_na)


def is_true(x: SpatRaster) -> SpatRaster:
    """Cell-wise ``isTRUE()``."""
    return _wrap(x, "is_true")


def is_false(x: SpatRaster) -> SpatRaster:
    """Cell-wise ``isFALSE()``."""
    return _wrap(x, "is_false")


def is_nan(x: SpatRaster) -> SpatRaster:
    """Cell-wise ``is.nan()``."""
    return _wrap(x, "isnan")


def is_finite(x: SpatRaster) -> SpatRaster:
    """Cell-wise ``is.finite()``."""
    return _wrap(x, "isfinite")


def is_infinite(x: SpatRaster) -> SpatRaster:
    """Cell-wise ``is.infinite()``."""
    return _wrap(x, "isinfinite")


def any_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """``anyNA()`` — any NA across layers."""
    return _wrap(x, "anynan", false_na)


def all_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """``allNA()`` — all values NA across layers."""
    return _wrap(x, "allnan", false_na)


def no_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """``noNA()`` — no NA values across layers."""
    return _wrap(x, "nonan", false_na)


def count_na(x: SpatRaster, n: int = 0) -> SpatRaster:
    """``countNA()`` — count NA values per cell across layers."""
    n = int(round(n))
    opt = _opt()
    if n == 1:
        r = x.is_wrapper("anynan", False, opt)
    else:
        r = x.countnan(n, opt)
    return messages(r, "countNA")


# ── Summary / reduction ──────────────────────────────────────────────────────

def _summarize(x: SpatRaster, fun: str, na_rm: bool = False,
               *extra: float, filename: str = "", **kw: Any) -> SpatRaster:
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    if extra:
        r = x.summary_numb(fun, list(extra), na_rm, opt)
    else:
        r = x.summary(fun, na_rm, opt)
    return messages(r, fun)


def which_max(x: SpatRaster) -> SpatRaster:
    """Layer index of maximum value — like R ``which.max()``."""
    return _summarize(x, "which.max", True)


def which_min(x: SpatRaster) -> SpatRaster:
    """Layer index of minimum value — like R ``which.min()``."""
    return _summarize(x, "which.min", True)


def which_lyr(x: SpatRaster) -> SpatRaster:
    """Which layer has a non-zero value — like R ``which.lyr()``."""
    return _summarize(x, "which", True)


def where_max(x: SpatRaster, values: bool = True) -> Any:
    """Cell and value of maximum — like R ``where.max()``."""
    out = x.where("max", values, _opt())
    messages(x, "where.max")
    return out


def where_min(x: SpatRaster, values: bool = True) -> Any:
    """Cell and value of minimum — like R ``where.min()``."""
    out = x.where("min", values, _opt())
    messages(x, "where.min")
    return out


def rast_sum(x: SpatRaster, *args: Any, na_rm: bool = False,
             filename: str = "", **kw: Any) -> SpatRaster:
    """Sum across layers — like R ``sum(SpatRaster)``."""
    return _summarize(x, "sum", na_rm, *args, filename=filename, **kw)


def rast_mean(x: SpatRaster, *args: Any, na_rm: bool = False,
              filename: str = "", **kw: Any) -> SpatRaster:
    """Mean across layers — like R ``mean(SpatRaster)``."""
    return _summarize(x, "mean", na_rm, *args, filename=filename, **kw)


def rast_min(x: SpatRaster, *args: Any, na_rm: bool = False,
             filename: str = "", **kw: Any) -> SpatRaster:
    """Min across layers — like R ``min(SpatRaster)``."""
    return _summarize(x, "min", na_rm, *args, filename=filename, **kw)


def rast_max(x: SpatRaster, *args: Any, na_rm: bool = False,
             filename: str = "", **kw: Any) -> SpatRaster:
    """Max across layers — like R ``max(SpatRaster)``."""
    return _summarize(x, "max", na_rm, *args, filename=filename, **kw)


def rast_median(x: SpatRaster, *args: Any, na_rm: bool = False,
                filename: str = "", **kw: Any) -> SpatRaster:
    """Median across layers — like R ``median(SpatRaster)``."""
    return _summarize(x, "median", na_rm, *args, filename=filename, **kw)


def stdev_rast(x: SpatRaster, *args: Any, pop: bool = True,
               na_rm: bool = False, filename: str = "", **kw: Any) -> SpatRaster:
    """Standard deviation across layers — like R ``stdev(SpatRaster)``."""
    fun = "std" if pop else "sd"
    return _summarize(x, fun, na_rm, *args, filename=filename, **kw)


def rast_modal(x: SpatRaster, *args: Any, ties: str = "first",
               na_rm: bool = False, filename: str = "", **kw: Any) -> SpatRaster:
    """Modal value across layers — like R ``modal(SpatRaster)``."""
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    extra = [float(a) for a in args]
    r = x.modal(extra, ties, na_rm, opt)
    return messages(r, "modal")


# ── Compare / logic functions ─────────────────────────────────────────────────

def compare_rast(
    x: SpatRaster,
    y: Union[SpatRaster, float, int],
    oper: str,
    false_na: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """
    Compare raster values — like R ``compare()``.

    *oper* is one of ``"=="  "!="  ">"  "<"  ">="  "<="``
    """
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    valid = {"==", "!=", ">", "<", ">=", "<="}
    if oper not in valid:
        raise ValueError(f"compare_rast: oper must be one of {valid}")
    if isinstance(y, SpatRaster):
        r = x.arith_rast(y, oper, false_na, opt)
    else:
        r = x.arith_numb([float(y)], oper, False, false_na, opt)
    return messages(r, oper)


def logic_rast_fn(
    x: SpatRaster,
    oper: str,
    false_na: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """
    Apply a logical mask operation — like R ``logic()``.

    *oper* is one of:
    ``"is.na"``  ``"not.na"``  ``"allNA"``  ``"anyNA"``  ``"noNA"``
    ``"is.infinite"``  ``"is.finite"``  ``"isTRUE"``  ``"isFALSE"``  ``"!"``
    """
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    mapping = {
        "is.na": "isnan",      "isNA": "isnan",
        "not.na": "isnotnan",  "not_na": "isnotnan",
        "allNA": "allnan",     "anyNA": "anynan",
        "noNA": "nonan",       "noneNA": "nonan",
        "is.infinite": "isinfinite",
        "is.finite": "isfinite",
        "isTRUE": "is_true",   "isFALSE": "is_false",
    }
    if oper == "!":
        r = x.arith_numb([0.0], "==", False, false_na, opt)
    elif oper in mapping:
        r = x.is_wrapper(mapping[oper], false_na, opt)
    else:
        raise ValueError(f"logic_rast_fn: unknown oper '{oper}'")
    return messages(r, "logic")


# ── Type coercion / type tests ────────────────────────────────────────────────

def as_int_rast(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Truncate to integer values — like R ``as.int(SpatRaster)``."""
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    r = x.math("trunc", opt)
    return messages(r, "as.int")


def as_bool_rast(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Convert to boolean (0/1) — like R ``as.bool(SpatRaster)``."""
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    r = x.is_wrapper("is_true", False, opt)
    return messages(r, "as.bool")


def is_bool_rast(x: SpatRaster) -> List[bool]:
    """True per layer where value type is boolean — like R ``is.bool()``."""
    return [v == 3 for v in x.valueType(False)]


def is_int_rast(x: SpatRaster) -> List[bool]:
    """True per layer where values are integers (non-categorical)."""
    vt = x.valueType(False)
    cats = x.hasCategories()
    return [(vt[i] == 1) and not cats[i] for i in range(len(vt))]


def is_num_rast(x: SpatRaster) -> List[bool]:
    """True per layer where values are numeric (non-categorical)."""
    vt = x.valueType(False)
    cats = x.hasCategories()
    return [(vt[i] < 2) and not cats[i] for i in range(len(vt))]


# ── Operator registration ─────────────────────────────────────────────────────

def register_operators() -> None:
    """
    Monkey-patch Python operators onto the pybind11 C++ types.

    Called once at import time by :mod:`terra.__init__`.
    """

    # ── SpatExtent ──
    SpatExtent.__add__  = _ext_add                             # type: ignore
    SpatExtent.__radd__ = lambda e, n: _ext_add(e, n)          # type: ignore
    SpatExtent.__sub__  = _ext_sub                             # type: ignore
    SpatExtent.__mul__  = _ext_mul                             # type: ignore
    SpatExtent.__rmul__ = lambda e, n: _ext_mul(e, n)          # type: ignore
    SpatExtent.__truediv__ = _ext_div                          # type: ignore
    SpatExtent.__mod__  = _ext_mod                             # type: ignore
    SpatExtent.__eq__   = lambda e, o: _ext_compare(e, o, "==")  # type: ignore
    SpatExtent.__ne__   = lambda e, o: _ext_compare(e, o, "!=")  # type: ignore
    SpatExtent.__lt__   = lambda e, o: _ext_compare(e, o, "<")   # type: ignore
    SpatExtent.__le__   = lambda e, o: _ext_compare(e, o, "<=")  # type: ignore
    SpatExtent.__gt__   = lambda e, o: _ext_compare(e, o, ">")   # type: ignore
    SpatExtent.__ge__   = lambda e, o: _ext_compare(e, o, ">=")  # type: ignore

    # ── SpatRaster ──
    for op in ("+", "-", "*", "/", "**", "//", "%"):
        _r_op = op

        def _make_op(o: str):
            def _op(a: SpatRaster, b: Any) -> SpatRaster:
                return _rast_binop(a, b, o)
            def _rop(a: SpatRaster, b: Any) -> SpatRaster:
                return _rast_rbinop(b, a, o)
            return _op, _rop

        fwd, rev = _make_op(_r_op)
        _sym_map = {
            "+": ("__add__", "__radd__"),
            "-": ("__sub__", "__rsub__"),
            "*": ("__mul__", "__rmul__"),
            "/": ("__truediv__", "__rtruediv__"),
            "**": ("__pow__", "__rpow__"),
            "//": ("__floordiv__", "__rfloordiv__"),
            "%": ("__mod__", "__rmod__"),
        }
        fa, ra = _sym_map[_r_op]
        setattr(SpatRaster, fa, fwd)
        setattr(SpatRaster, ra, rev)

    # Unary -r  →  0 - r
    SpatRaster.__neg__ = lambda r: _rast_rbinop(0.0, r, "-")  # type: ignore
    # ~r  →  r == 0
    SpatRaster.__invert__ = lambda r: _rast_binop(r, 0.0, "==")  # type: ignore

    # Comparison
    for op in ("==", "!=", "<", "<=", ">", ">="):
        _c_op = op

        def _make_cmp(o: str):
            def _cmp(a: SpatRaster, b: Any) -> SpatRaster:
                return _rast_binop(a, b, o)
            return _cmp

        _sym_cmp = {
            "==": "__eq__", "!=": "__ne__",
            "<":  "__lt__", "<=": "__le__",
            ">":  "__gt__", ">=": "__ge__",
        }
        setattr(SpatRaster, _sym_cmp[_c_op], _make_cmp(_c_op))

    # Logic
    SpatRaster.__and__ = lambda a, b: _rast_logic_op(a, b, "&")   # type: ignore
    SpatRaster.__or__  = lambda a, b: _rast_logic_op(a, b, "|")   # type: ignore
    SpatRaster.__rand__ = lambda a, b: _rast_logic_op(a, bool(b), "&")  # type: ignore
    SpatRaster.__ror__  = lambda a, b: _rast_logic_op(a, bool(b), "|")  # type: ignore

    # ── SpatVector ──
    SpatVector.__add__ = lambda a, b: _vect_binop(a, b, "+")  # type: ignore
    SpatVector.__mul__ = lambda a, b: _vect_binop(a, b, "*")  # type: ignore
    SpatVector.__sub__ = lambda a, b: _vect_binop(a, b, "-")  # type: ignore
