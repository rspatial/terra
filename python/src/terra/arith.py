"""
Arithmetic, comparison, logical, and NA/summary operators for terra objects.

Operator overloads are registered at import time via :func:`register_operators`.

**SpatRaster** supports: ``+  -  *  /  //  %  **`` (cell-wise, with another
SpatRaster, a scalar, a list of scalars, or a 2-D NumPy array);
comparison ``== != < <= > >=``; logical ``& |``; and unary ``-`` and ``~``.

**SpatExtent** supports: ``+ - * / %`` with a numeric scalar (expand, shrink,
scale from centre, and align), ``+`` and ``*`` with another SpatExtent (union
and intersection), and ``/ SpatExtent`` (returns the width/height ratio);
comparison ``== != < <= > >=``.

**SpatVector** supports: ``+ * -`` between two SpatVectors (union, intersect,
erase).
"""

from __future__ import annotations

from typing import Any, List, Union

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


# ── Internal helpers ─────────────────────────────────────────────────────────

def _opt() -> SpatOptions:
    return SpatOptions()


def _to_floats(x: Any) -> List[float]:
    if isinstance(x, (int, float)):
        return [float(x)]
    return [float(v) for v in x]


# ── SpatExtent arithmetic ────────────────────────────────────────────────────

def _ext_add(e: SpatExtent, n: Any) -> SpatExtent:
    from .extent import ext as _ext
    if isinstance(n, SpatExtent):
        ec = e.deepcopy()
        ec.union(n)
        return ec
    v = e.vector
    ns = _to_floats(n)
    while len(ns) < 4:
        ns.extend(ns)
    ns = ns[:4]
    return _ext(v[0] - ns[0], v[1] + ns[1], v[2] - ns[2], v[3] + ns[3])


def _ext_sub(e: SpatExtent, n: Any) -> SpatExtent:
    from .extent import ext as _ext
    v = e.vector
    ns = _to_floats(n)
    while len(ns) < 4:
        ns.extend(ns)
    ns = ns[:4]
    result = _ext(v[0] + ns[0], v[1] - ns[1], v[2] + ns[2], v[3] - ns[3])
    if not result.valid:
        raise ValueError("- : this would create an invalid extent")
    return result


def _ext_mul(e: SpatExtent, n: Any) -> Any:
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


# ── NA / logical cell-wise tests ─────────────────────────────────────────────

def _wrap(r: SpatRaster, name: str, false_na: bool = False) -> SpatRaster:
    r = r.is_wrapper(name, false_na, _opt())
    return messages(r, name)


def is_na(x: SpatRaster) -> SpatRaster:
    """
    Create a boolean raster indicating which cells are NA.

    Computations are cell-wise. The terra package does not distinguish between
    NA (not available) and NaN (not a number); in most cases this state is
    represented by NaN.

    Args:
        x: SpatRaster to test.

    Returns:
        SpatRaster with TRUE (1) where ``x`` is NA and FALSE (0) elsewhere.
    """
    return _wrap(x, "isnan")


def not_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """
    Create a boolean raster indicating which cells are not NA.

    Shortcut to avoid the two-step ``~is_na(x)``.

    Args:
        x: SpatRaster to test.
        false_na: If True, cells that would otherwise be FALSE are set to NA
            instead.

    Returns:
        SpatRaster with TRUE where ``x`` is not NA.
    """
    return _wrap(x, "isnotnan", false_na)


def is_true(x: SpatRaster) -> SpatRaster:
    """
    Create a boolean raster that is TRUE for cells with a non-zero, non-NA value.

    Equivalent to ``as_bool_rast(x)``.

    Args:
        x: SpatRaster to test.

    Returns:
        SpatRaster with boolean values.
    """
    return _wrap(x, "is_true")


def is_false(x: SpatRaster) -> SpatRaster:
    """
    Create a boolean raster that is TRUE for cells with the value zero.

    Equivalent to ``~as_bool_rast(x)``.

    Args:
        x: SpatRaster to test.

    Returns:
        SpatRaster with boolean values.
    """
    return _wrap(x, "is_false")


def is_nan(x: SpatRaster) -> SpatRaster:
    """
    Create a boolean raster indicating which cells are NaN/NA.

    The terra package does not distinguish between NA and NaN; this method
    behaves identically to :func:`is_na`.

    Args:
        x: SpatRaster to test.

    Returns:
        SpatRaster with TRUE (1) where ``x`` is NaN/NA.
    """
    return _wrap(x, "isnan")


def is_finite(x: SpatRaster) -> SpatRaster:
    """
    Create a boolean raster indicating which cells are finite.

    Args:
        x: SpatRaster to test.

    Returns:
        SpatRaster with TRUE (1) where ``x`` is finite.
    """
    return _wrap(x, "isfinite")


def is_infinite(x: SpatRaster) -> SpatRaster:
    """
    Create a boolean raster indicating which cells are infinite.

    Args:
        x: SpatRaster to test.

    Returns:
        SpatRaster with TRUE (1) where ``x`` is infinite.
    """
    return _wrap(x, "isinfinite")


def any_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """
    Create a boolean raster that is TRUE where any layer has an NA value.

    Computations are across layers; each output cell is TRUE if at least one
    layer value is NA.

    Args:
        x: SpatRaster.
        false_na: If True, cells that would otherwise be FALSE are set to NA.

    Returns:
        SpatRaster (single layer).
    """
    return _wrap(x, "anynan", false_na)


def all_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """
    Create a boolean raster that is TRUE where all layers have NA values.

    Args:
        x: SpatRaster.
        false_na: If True, cells that would otherwise be FALSE are set to NA.

    Returns:
        SpatRaster (single layer).
    """
    return _wrap(x, "allnan", false_na)


def no_na(x: SpatRaster, false_na: bool = False) -> SpatRaster:
    """
    Create a boolean raster that is TRUE where no layer has an NA value.

    Args:
        x: SpatRaster.
        false_na: If True, cells that would otherwise be FALSE are set to NA.

    Returns:
        SpatRaster (single layer).
    """
    return _wrap(x, "nonan", false_na)


def count_na(x: SpatRaster, n: int = 0) -> SpatRaster:
    """
    Count the number of NA values across layers for each cell.

    Args:
        x: SpatRaster.
        n: If ``n > 0``, cell values are TRUE if at least ``n`` of the layers
            are NA, rather than returning a count.

    Returns:
        SpatRaster with the NA count (or a boolean indicator when ``n > 0``).
    """
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
    """
    Return the layer index of the maximum value for each cell.

    Args:
        x: SpatRaster with one or more layers.

    Returns:
        Single-layer SpatRaster with the 1-based layer index of the maximum.
    """
    return _summarize(x, "which.max", True)


def which_min(x: SpatRaster) -> SpatRaster:
    """
    Return the layer index of the minimum value for each cell.

    Args:
        x: SpatRaster with one or more layers.

    Returns:
        Single-layer SpatRaster with the 1-based layer index of the minimum.
    """
    return _summarize(x, "which.min", True)


def which_lyr(x: SpatRaster) -> SpatRaster:
    """
    Return the layer index of the first non-zero value for each cell.

    Args:
        x: SpatRaster with one or more layers.

    Returns:
        Single-layer SpatRaster.
    """
    return _summarize(x, "which", True)


def where_max(x: SpatRaster, values: bool = True) -> Any:
    """
    Return the cell numbers for cells with the maximum value per layer.

    Args:
        x: SpatRaster.
        values: If True, the maximum cell values are also included in the
            output alongside the cell numbers.

    Returns:
        A list (one element per layer) of arrays containing cell numbers and,
        if ``values=True``, the maximum values.
    """
    out = x.where("max", values, _opt())
    messages(x, "where.max")
    return out


def where_min(x: SpatRaster, values: bool = True) -> Any:
    """
    Return the cell numbers for cells with the minimum value per layer.

    Args:
        x: SpatRaster.
        values: If True, the minimum cell values are also included in the
            output alongside the cell numbers.

    Returns:
        A list (one element per layer) of arrays containing cell numbers and,
        if ``values=True``, the minimum values.
    """
    out = x.where("min", values, _opt())
    messages(x, "where.min")
    return out


def rast_sum(x: SpatRaster, *args: Any, na_rm: bool = False,
             filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the sum across layers for each cell.

    Additional SpatRasters or numeric values passed as positional arguments
    are included in the computation. Each numeric value is treated as a
    constant-valued layer.

    Args:
        x: SpatRaster.
        *args: Additional SpatRasters or numeric values.
        na_rm: If True, NA values are ignored. If False, NA is returned
            whenever any layer value is NA.
        filename: Output filename (empty string for in-memory result).
        **kw: Additional arguments passed to the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    return _summarize(x, "sum", na_rm, *args, filename=filename, **kw)


def rast_mean(x: SpatRaster, *args: Any, na_rm: bool = False,
              filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the mean across layers for each cell.

    Args:
        x: SpatRaster.
        *args: Additional SpatRasters or numeric values included in the mean.
        na_rm: If True, NA values are ignored.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    return _summarize(x, "mean", na_rm, *args, filename=filename, **kw)


def rast_min(x: SpatRaster, *args: Any, na_rm: bool = False,
             filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the minimum across layers for each cell.

    Args:
        x: SpatRaster.
        *args: Additional SpatRasters or numeric values.
        na_rm: If True, NA values are ignored.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    return _summarize(x, "min", na_rm, *args, filename=filename, **kw)


def rast_max(x: SpatRaster, *args: Any, na_rm: bool = False,
             filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the maximum across layers for each cell.

    Args:
        x: SpatRaster.
        *args: Additional SpatRasters or numeric values.
        na_rm: If True, NA values are ignored.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    return _summarize(x, "max", na_rm, *args, filename=filename, **kw)


def rast_median(x: SpatRaster, *args: Any, na_rm: bool = False,
                filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the median across layers for each cell.

    Args:
        x: SpatRaster.
        *args: Additional numeric values treated as constant layers.
        na_rm: If True, NA values are ignored. Must be a boolean value.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    return _summarize(x, "median", na_rm, *args, filename=filename, **kw)


def stdev_rast(x: SpatRaster, *args: Any, pop: bool = True,
               na_rm: bool = False, filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the standard deviation across layers for each cell.

    If ``pop=True`` the population standard deviation is computed:
    ``sqrt(sum((x - mean(x))^2) / n)``. This differs from the sample standard
    deviation (which uses ``n - 1`` as denominator).

    Args:
        x: SpatRaster.
        *args: Additional SpatRasters or numeric values included in the
            computation.
        pop: If True (default), compute the population standard deviation.
            If False, compute the sample standard deviation.
        na_rm: If True, NA values are ignored.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    fun = "std" if pop else "sd"
    return _summarize(x, fun, na_rm, *args, filename=filename, **kw)


def rast_modal(x: SpatRaster, *args: Any, ties: str = "first",
               na_rm: bool = False, filename: str = "", **kw: Any) -> SpatRaster:
    """
    Compute the mode (most frequent value) across layers for each cell.

    Args:
        x: SpatRaster.
        *args: Additional SpatRasters or numeric values to include.
        ties: How to handle ties among equally frequent values. One of
            ``"first"``, ``"last"``, ``"random"``, ``"lowest"``,
            ``"highest"``, or ``"NA"``.
        na_rm: If True, NA values are ignored.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        Single-layer SpatRaster.
    """
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    extra = [float(a) for a in args]
    r = x.modal(extra, ties, na_rm, opt)
    return messages(r, "modal")


# ── Compare / logic ───────────────────────────────────────────────────────────

def compare_rast(
    x: SpatRaster,
    y: Union[SpatRaster, float, int],
    oper: str,
    false_na: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """
    Compare cell values using a relational operator.

    Computations are cell-wise. If multiple SpatRasters are used they must
    have the same geometry (extent and resolution).

    Unlike the ``==`` operator, this function can return NA instead of FALSE
    via ``false_na``, and supports writing to a file.

    Args:
        x: SpatRaster.
        y: SpatRaster or numeric value to compare against.
        oper: Comparison operator string: one of ``"=="``, ``"!="``, ``">"``,
            ``"<"``, ``">="``, ``"<="``.
        false_na: If True, cells that compare as FALSE become NA instead.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        SpatRaster with boolean (0/1/NA) values.
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
    Apply a logical mask or NA-detection operation to a raster.

    Unlike the operator forms, this function can return NA instead of FALSE
    via ``false_na``, and supports writing to a file.

    Args:
        x: SpatRaster.
        oper: Operation name. One of:

            * ``"!"``         — logical NOT (equivalent to ``x == 0``)
            * ``"is.na"``     — TRUE where cell is NA
            * ``"not.na"``    — TRUE where cell is not NA
            * ``"allNA"``     — TRUE where all layers are NA
            * ``"anyNA"``     — TRUE where any layer is NA
            * ``"noNA"`` / ``"noneNA"`` — TRUE where no layer is NA
            * ``"is.infinite"``
            * ``"is.finite"``
            * ``"isTRUE"``    — TRUE where cell is non-zero and non-NA
            * ``"isFALSE"``   — TRUE where cell is zero

        false_na: If True, cells that would be FALSE become NA instead.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        SpatRaster with boolean (0/1/NA) values.
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


# ── Type coercion / type queries ──────────────────────────────────────────────

def as_int_rast(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """
    Force raster values to integer by truncation.

    In-memory values remain stored as numeric but the layer is flagged as
    integer and written as an integer data type when saved to file (if the
    format supports it).

    Args:
        x: SpatRaster.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        SpatRaster with integer-typed values.
    """
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    r = x.math("trunc", opt)
    return messages(r, "as.int")


def as_bool_rast(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """
    Force raster values to boolean (TRUE/FALSE).

    Cells with a non-zero, non-NA value become TRUE (1); zero cells become
    FALSE (0); NA cells remain NA.

    Args:
        x: SpatRaster.
        filename: Output filename.
        **kw: Additional arguments for the file writer.

    Returns:
        SpatRaster with boolean-typed values.
    """
    from .generics import _opt as _g_opt
    opt = _g_opt(filename, **kw)
    r = x.is_wrapper("is_true", False, opt)
    return messages(r, "as.bool")


def is_bool_rast(x: SpatRaster) -> List[bool]:
    """
    Test whether each layer holds boolean values.

    Args:
        x: SpatRaster.

    Returns:
        List of bool, one per layer — True if the layer is boolean-typed.
    """
    return [v == 3 for v in x.valueType(False)]


def is_int_rast(x: SpatRaster) -> List[bool]:
    """
    Test whether each layer holds integer values (non-categorical).

    Args:
        x: SpatRaster.

    Returns:
        List of bool, one per layer — True if the layer is integer-typed and
        not categorical.
    """
    vt = x.valueType(False)
    cats = x.hasCategories()
    return [(vt[i] == 1) and not cats[i] for i in range(len(vt))]


def is_num_rast(x: SpatRaster) -> List[bool]:
    """
    Test whether each layer holds numeric (floating-point) values.

    Args:
        x: SpatRaster.

    Returns:
        List of bool, one per layer — True if the layer is numeric and not
        categorical.
    """
    vt = x.valueType(False)
    cats = x.hasCategories()
    return [(vt[i] < 2) and not cats[i] for i in range(len(vt))]


# ── Operator registration ─────────────────────────────────────────────────────

def register_operators() -> None:
    """
    Attach Python arithmetic, comparison, and logical operators to the C++ types.

    Called once at import time by ``terra/__init__.py``. After this call,
    standard Python operators (``+``, ``-``, ``*``, ``/``, ``==``, ``&``,
    etc.) work directly on SpatRaster, SpatVector, and SpatExtent instances.
    """

    # ── SpatExtent ──
    SpatExtent.__add__     = _ext_add                               # type: ignore
    SpatExtent.__radd__    = lambda e, n: _ext_add(e, n)            # type: ignore
    SpatExtent.__sub__     = _ext_sub                               # type: ignore
    SpatExtent.__mul__     = _ext_mul                               # type: ignore
    SpatExtent.__rmul__    = lambda e, n: _ext_mul(e, n)            # type: ignore
    SpatExtent.__truediv__ = _ext_div                               # type: ignore
    SpatExtent.__mod__     = _ext_mod                               # type: ignore
    SpatExtent.__eq__      = lambda e, o: _ext_compare(e, o, "==") # type: ignore
    SpatExtent.__ne__      = lambda e, o: _ext_compare(e, o, "!=") # type: ignore
    SpatExtent.__lt__      = lambda e, o: _ext_compare(e, o, "<")  # type: ignore
    SpatExtent.__le__      = lambda e, o: _ext_compare(e, o, "<=") # type: ignore
    SpatExtent.__gt__      = lambda e, o: _ext_compare(e, o, ">")  # type: ignore
    SpatExtent.__ge__      = lambda e, o: _ext_compare(e, o, ">=") # type: ignore

    # ── SpatRaster ──
    _sym_map = {
        "+":  ("__add__",       "__radd__"),
        "-":  ("__sub__",       "__rsub__"),
        "*":  ("__mul__",       "__rmul__"),
        "/":  ("__truediv__",   "__rtruediv__"),
        "**": ("__pow__",       "__rpow__"),
        "//": ("__floordiv__",  "__rfloordiv__"),
        "%":  ("__mod__",       "__rmod__"),
    }
    for op, (fa, ra) in _sym_map.items():
        def _make_op(o: str):
            def _fwd(a: SpatRaster, b: Any) -> SpatRaster:
                return _rast_binop(a, b, o)
            def _rev(a: SpatRaster, b: Any) -> SpatRaster:
                return _rast_rbinop(b, a, o)
            return _fwd, _rev
        fwd, rev = _make_op(op)
        setattr(SpatRaster, fa, fwd)
        setattr(SpatRaster, ra, rev)

    SpatRaster.__neg__    = lambda r: _rast_rbinop(0.0, r, "-")    # type: ignore
    SpatRaster.__invert__ = lambda r: _rast_binop(r, 0.0, "==")    # type: ignore

    _cmp_map = {
        "==": "__eq__", "!=": "__ne__",
        "<":  "__lt__", "<=": "__le__",
        ">":  "__gt__", ">=": "__ge__",
    }
    for op, dunder in _cmp_map.items():
        def _make_cmp(o: str):
            def _cmp(a: SpatRaster, b: Any) -> SpatRaster:
                return _rast_binop(a, b, o)
            return _cmp
        setattr(SpatRaster, dunder, _make_cmp(op))

    SpatRaster.__and__  = lambda a, b: _rast_logic_op(a, b, "&")   # type: ignore
    SpatRaster.__or__   = lambda a, b: _rast_logic_op(a, b, "|")   # type: ignore
    SpatRaster.__rand__ = lambda a, b: _rast_logic_op(a, bool(b), "&")  # type: ignore
    SpatRaster.__ror__  = lambda a, b: _rast_logic_op(a, bool(b), "|")  # type: ignore

    # ── SpatVector ──
    SpatVector.__add__ = lambda a, b: _vect_binop(a, b, "+")       # type: ignore
    SpatVector.__mul__ = lambda a, b: _vect_binop(a, b, "*")       # type: ignore
    SpatVector.__sub__ = lambda a, b: _vect_binop(a, b, "-")       # type: ignore
