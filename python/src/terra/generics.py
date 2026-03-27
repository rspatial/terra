"""
High-level functions mirroring R ``terra`` generics.

Each function corresponds closely to the same-named R generic so that R
workflows translate with minimal renaming.  The underlying C++ methods are
called directly on the pybind11-bound objects; no R S4 dispatch needed.
"""

from __future__ import annotations

from typing import Any, List, Optional, Sequence, Union

from ._helpers import messages
from ._terra import SpatExtent, SpatOptions, SpatRaster, SpatVector

# Capture raw C++ method references at import time, before register_methods()
# patches them.  Every Python-level wrapper that delegates to a C++ method
# with the same name must use these saved references to avoid infinite
# recursion when the method-style API (r.crop(e), etc.) is in use.
_cpp = {
    "rast.crop":       SpatRaster.crop,
    "rast.classify":   SpatRaster.classify,
    "rast.boundaries": SpatRaster.boundaries,
    "rast.patches":    SpatRaster.patches,
    "rast.terrain":    SpatRaster.terrain,
    "rast.sieve":      SpatRaster.sieve,
    "rast.stretch":    SpatRaster.stretch,
    "rast.trim":       SpatRaster.trim,
    "rast.flip":       SpatRaster.flip,
    "rast.rotate":     SpatRaster.rotate,
    "rast.shift":      SpatRaster.shift,
}

__all__ = [
    # dimensions / metadata
    "nrow", "ncol", "nlyr", "ncell", "res", "origin",
    # helpers
    "spat_options", "deepcopy", "tighten",
    # extent
    "ext_align",
    # raster geometry
    "is_rotated", "is_flipped", "flip", "rotate", "shift", "rescale",
    "trans", "trim", "rev_raster",
    # raster values
    "clamp", "clamp_ts", "classify", "subst", "cover", "diff_raster",
    "disagg", "segregate", "selectRange", "sort_raster",
    "range_fill", "weighted_mean",
    # raster analysis
    "boundaries", "patches", "cellSize", "surfArea", "terrain",
    "sieve", "rectify", "stretch", "scale_linear", "scale_raster",
    "quantile_raster", "atan_2",
    # raster processing
    "crop", "mask", "project_raster", "resample", "intersect_rast",
    # vector
    "project_vector", "shift_vect", "rotate_vect", "rescale_vect",
    "trans_vect",
    # scoff
    "scoff", "scoff_set",
]


def _opt(filename: str = "", **kw: Any) -> SpatOptions:
    """Create a SpatOptions, applying common keyword overrides."""
    opt = SpatOptions()
    if filename:
        opt.filenames = [filename]
    for k, v in kw.items():
        setattr(opt, k, v)
    return opt


# ── Dimensions ───────────────────────────────────────────────────────────────

def nrow(x: Any) -> int:
    """Number of rows."""
    return x.nrow()


def ncol(x: Any) -> int:
    """Number of columns."""
    return x.ncol()


def nlyr(x: Any) -> int:
    """Number of layers."""
    return x.nlyr()


def ncell(x: Any) -> int:
    """Total number of cells."""
    return x.ncell()


def res(x: SpatRaster) -> List[float]:
    """Spatial resolution (x, y)."""
    return x.res()


def origin(x: SpatRaster) -> List[float]:
    """Cell origin (xmin % xres, ymin % yres)."""
    return x.origin


# ── Options / utilities ──────────────────────────────────────────────────────

def spat_options() -> SpatOptions:
    """Like R ``spatOptions()`` — return a new default :class:`SpatOptions`."""
    return SpatOptions()


def deepcopy(x: Any) -> Any:
    """Deep copy — like R ``deepcopy()``."""
    return x.deepcopy()


def tighten(x: SpatRaster) -> SpatRaster:
    """Collapse multiple sources — like R ``tighten()``."""
    x = x.collapse_sources()
    return messages(x, "tighten")


# ── Extent align ─────────────────────────────────────────────────────────────

def ext_align(e: SpatExtent, y: Any, snap: str = "near") -> SpatExtent:
    """
    Align a :class:`SpatExtent` to raster *y*, like R ``align(x, y)``.

    *y* can be a :class:`SpatRaster` (use ``snap``) or a ``float`` resolution.
    """
    if isinstance(y, SpatRaster):
        return y.align(e, snap.lower())
    # numeric resolution
    return e.align(float(y), "")


# ── Raster geometry ──────────────────────────────────────────────────────────

def is_rotated(x: SpatRaster) -> bool:
    """True if the raster has a rotation angle."""
    return x.is_rotated()


def is_flipped(x: SpatRaster) -> bool:
    """True if the raster is flipped (top-to-bottom)."""
    return x.is_flipped()


def flip(x: SpatRaster, direction: str = "vertical", filename: str = "", **kw: Any) -> SpatRaster:
    """
    Flip a raster — like R ``flip()``.

    *direction*: ``"vertical"`` (default) or ``"horizontal"``.
    """
    d = direction.lower()
    if d not in ("vertical", "horizontal"):
        raise ValueError("direction must be 'vertical' or 'horizontal'")
    opt = _opt(filename, **kw)
    x = _cpp["rast.flip"](x, d == "vertical", opt)
    return messages(x, "flip")


def rotate(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Rotate a raster 180° (i.e. shift longitude 0/360) — like R ``rotate()``."""
    opt = _opt(filename, **kw)
    x = _cpp["rast.rotate"](x, True, opt)
    return messages(x, "rotate")


def shift(
    x: Any,
    dx: float = 0.0,
    dy: float = 0.0,
    filename: str = "",
    **kw: Any,
) -> Any:
    """
    Shift coordinates — like R ``shift()``.

    Works for :class:`SpatRaster`, :class:`SpatExtent`, or :class:`SpatVector`.
    """
    if isinstance(x, SpatRaster):
        opt = _opt(filename, **kw)
        x = _cpp["rast.shift"](x, dx, dy, opt)
        return messages(x, "shift")
    if isinstance(x, SpatExtent):
        v = x.vector
        from .extent import ext as _ext
        return _ext(v[0] + dx, v[1] + dx, v[2] + dy, v[3] + dy)
    if isinstance(x, SpatVector):
        x = x.shift(dx, dy)
        return messages(x, "shift")
    raise TypeError("shift: unsupported type")


def rescale(
    x: Any,
    fx: float = 0.5,
    fy: Optional[float] = None,
    x0: Optional[float] = None,
    y0: Optional[float] = None,
) -> Any:
    """
    Rescale — like R ``rescale()``.

    Works for :class:`SpatRaster` and :class:`SpatVector`.
    """
    if fy is None:
        fy = fx
    if isinstance(x, SpatRaster):
        e = x.extent.vector
        _x0 = x0 if x0 is not None else (e[0] + e[1]) / 2.0
        _y0 = y0 if y0 is not None else (e[2] + e[3]) / 2.0
        ex = [_x0 + fx * (e[0] - _x0), _x0 + fx * (e[1] - _x0)]
        ey = [_y0 + fy * (e[2] - _y0), _y0 + fy * (e[3] - _y0)]
        from ._terra import SpatExtent
        new_ext = SpatExtent(ex[0], ex[1], ey[0], ey[1])
        r = x.deepcopy()
        r.extent = new_ext
        return messages(r, "rescale")
    if isinstance(x, SpatVector):
        e = x.extent().vector
        _x0 = x0 if x0 is not None else (e[0] + e[1]) / 2.0
        _y0 = y0 if y0 is not None else (e[2] + e[3]) / 2.0
        x = x.rescale(fx, fy, float(_x0), float(_y0))
        return messages(x, "rescale")
    raise TypeError("rescale: unsupported type")


def trans(x: Any, filename: str = "", **kw: Any) -> Any:
    """Transpose rows/columns — like R ``trans()`` / ``t()``."""
    if isinstance(x, SpatRaster):
        opt = _opt(filename, **kw)
        x = x.transpose(opt)
        return messages(x, "trans")
    if isinstance(x, SpatVector):
        x = x.transpose()
        return messages(x, "trans")
    raise TypeError("trans: unsupported type")


def trim(
    x: SpatRaster,
    padding: int = 0,
    value: float = float("nan"),
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Trim NA margins — like R ``trim()``."""
    import math

    opt = _opt(filename, **kw)
    na = float("nan") if math.isnan(value) else float(value)
    x = _cpp["rast.trim"](x, na, int(padding), opt)
    return messages(x, "trim")


def rev_raster(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Reverse layer order — like R ``rev()`` on a SpatRaster."""
    opt = _opt(filename, **kw)
    x = x.reverse(opt)
    return messages(x, "rev")


# ── Raster values ─────────────────────────────────────────────────────────────

def range_fill(
    x: SpatRaster,
    limit: int,
    circular: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Fill range — like R ``rangeFill()``."""
    opt = _opt(filename, **kw)
    x = x.fill_range(int(limit), circular, opt)
    return messages(x, "rangeFill")


def weighted_mean(
    x: SpatRaster,
    w: Union[SpatRaster, List[float]],
    na_rm: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Weighted mean — like R ``weighted.mean()``."""
    opt = _opt(filename, **kw)
    if isinstance(w, SpatRaster):
        x = x.wmean_rast(w, na_rm, opt)
    else:
        x = x.wmean_vect(list(w), na_rm, opt)
    return messages(x, "weighted.mean")


def clamp(
    x: SpatRaster,
    lower: Union[float, SpatRaster] = float("-inf"),
    upper: Union[float, SpatRaster] = float("inf"),
    values: bool = True,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Clamp values to [lower, upper] — like R ``clamp()``."""
    opt = _opt(filename, **kw)
    empty = SpatRaster()
    rlow = isinstance(lower, SpatRaster)
    rupp = isinstance(upper, SpatRaster)
    if rlow and rupp:
        x = x.clamp_raster(lower, upper, [float("nan")], [float("nan")], values, opt)
    elif rlow:
        x = x.clamp_raster(lower, empty, [float("nan")], [float(upper)], values, opt)
    elif rupp:
        x = x.clamp_raster(empty, upper, [float(lower)], [float("nan")], values, opt)
    else:
        x = x.clamp_raster(empty, empty, [float(lower)], [float(upper)], values, opt)
    return messages(x, "clamp")


def clamp_ts(
    x: SpatRaster,
    min: bool = False,
    max: bool = True,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Clamp time series to min/max — like R ``clamp_ts()``."""
    opt = _opt(filename, **kw)
    x = x.clamp_ts(min, max, opt)
    return messages(x, "clamp_ts")


def classify(
    x: SpatRaster,
    rcl: Any,
    include_lowest: bool = False,
    right: Optional[bool] = True,
    others: Optional[float] = None,
    brackets: bool = True,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Reclassify values — like R ``classify()``."""
    opt = _opt(filename, **kw)

    try:
        import numpy as np  # type: ignore

        if isinstance(rcl, np.ndarray):
            rcl = rcl.tolist()
    except ImportError:
        pass

    rcl = [list(row) for row in rcl]
    ncols = len(rcl[0]) if rcl else 0

    if ncols == 2:
        from_v = [float(r[0]) for r in rcl]
        to_v = [float(r[1]) for r in rcl]
        use_others = others is not None
        ov = float(others) if use_others else float("nan")
        x = x.lookup_classify(from_v, to_v, use_others, ov, opt)
        return messages(x, "classify")

    right_i = 2 if right is None else (1 if right else 0)
    il = bool(include_lowest)
    use_others = others is not None
    ov = float(others) if use_others else 0.0
    flat = [float(v) for row in rcl for v in row]
    x = _cpp["rast.classify"](x, flat, ncols, right_i, il, use_others, ov, False, brackets, False, opt)
    return messages(x, "classify")


def subst(
    x: SpatRaster,
    from_val: Any,
    to_val: Any,
    others: Optional[float] = None,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Substitute values — like R ``subst()``."""
    opt = _opt(filename, **kw)
    use_others = others is not None
    ov = float(others) if use_others else float("nan")
    from_v = [float(v) for v in from_val]
    to_v = [float(v) for v in to_val]
    x = x.lookup_subst(from_v, to_v, use_others, ov, opt)
    return messages(x, "subst")


def cover(
    x: SpatRaster,
    y: Optional[SpatRaster] = None,
    values: Optional[float] = None,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Replace NA with values from *y* — like R ``cover()``."""
    opt = _opt(filename, **kw)
    na_vals = [float("nan")] if values is None else [float(values)]
    if y is None:
        x = x.cover_self(na_vals, opt)
    else:
        x = x.cover(y, na_vals, opt)
    return messages(x, "cover")


def diff_raster(x: SpatRaster, lag: int = 1, filename: str = "", **kw: Any) -> SpatRaster:
    """Layer-to-layer differences — like R ``diff()``."""
    n = x.nlyr()
    lag = int(round(lag))
    if lag < 1 or lag >= n:
        raise ValueError("diff: lag must be > 0 and < nlyr(x)")
    opt = _opt(filename, **kw)
    # y = first (n-lag) layers, x = last (n-lag) layers
    keep_y = list(range(n - lag))
    keep_x = list(range(lag, n))
    y_r = x.subset(keep_y, SpatOptions())
    x_r = x.subset(keep_x, SpatOptions())
    x_r = x_r.arith_rast(y_r, "-", False, opt)
    return messages(x_r, "diff")


def disagg(
    x: SpatRaster,
    fact: Union[int, Sequence[int]],
    method: str = "near",
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Disaggregate (increase resolution) — like R ``disagg()``."""
    method = method.lower()
    if method not in ("near", "bilinear"):
        raise ValueError("disagg: method must be 'near' or 'bilinear'")
    if method == "bilinear":
        # disagg geometry then resample
        if isinstance(fact, int):
            f = [fact, fact]
        else:
            f = list(fact)[:2]
        template = x.set_resolution(x.res()[0] / f[0], x.res()[1] / f[1])
        return resample(x, template, method="bilinear", filename=filename, **kw)
    opt = _opt(filename, **kw)
    if isinstance(fact, int):
        fact_list = [fact, fact]
    else:
        fact_list = [int(f) for f in list(fact)[:2]]
    x = x.disaggregate(fact_list, opt)
    return messages(x, "disagg")


def segregate(
    x: SpatRaster,
    classes: Optional[List[float]] = None,
    keep: bool = False,
    other: float = 0.0,
    round_vals: bool = False,
    digits: int = 0,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Segregate values into binary layers — like R ``segregate()``."""
    opt = _opt(filename, **kw)
    cls = classes if classes is not None else []
    x = x.separate(cls, keep, float(other), round_vals, digits, opt)
    return messages(x, "segregate")


def selectRange(
    x: SpatRaster,
    y: SpatRaster,
    z: int = 1,
    repint: int = 0,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Select values from a range of layers — like R ``selectRange()``."""
    opt = _opt(filename, **kw)
    x = x.selRange(y, z, repint, opt)
    return messages(x, "selectRange")


def sort_raster(
    x: SpatRaster,
    decreasing: bool = False,
    order: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Sort layer values — like R ``sort()`` on SpatRaster."""
    opt = _opt(filename, **kw)
    x = x.sort(decreasing, order, opt)
    return messages(x, "sort")


# ── Raster analysis ──────────────────────────────────────────────────────────

def boundaries(
    x: SpatRaster,
    classes: bool = False,
    inner: bool = True,
    directions: int = 8,
    falseval: float = 0.0,
    ignore_na: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Edge detection — like R ``boundaries()``."""
    opt = _opt(filename, **kw)
    btype = "inner" if inner else "outer"
    x = _cpp["rast.boundaries"](x, classes, ignore_na, btype, directions, falseval, opt)
    return messages(x, "boundaries")


def patches(
    x: SpatRaster,
    directions: int = 4,
    zero_as_na: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Label connected regions — like R ``patches()``."""
    opt = _opt(filename, **kw)
    x = _cpp["rast.patches"](x, directions, zero_as_na, opt)
    return messages(x, "patches")


def cellSize(
    x: SpatRaster,
    mask: bool = False,
    unit: str = "m",
    transform: bool = True,
    rcx: int = 100,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Cell area — like R ``cellSize()``."""
    opt = _opt(filename, **kw)
    x = x.rst_area(mask, unit, transform, rcx, opt)
    return messages(x, "cellSize")


def surfArea(x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Surface area — like R ``surfArea()``."""
    opt = _opt(filename, **kw)
    x = x.surface_area(opt)
    return messages(x, "surfArea")


def terrain(
    x: SpatRaster,
    v: Union[str, List[str]] = "slope",
    neighbors: int = 8,
    unit: str = "degrees",
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Terrain attributes — like R ``terrain()``."""
    unit = unit.lower()
    if unit not in ("degrees", "radians"):
        raise ValueError("unit must be 'degrees' or 'radians'")
    opt = _opt(filename, **kw)
    vlist = [v] if isinstance(v, str) else list(v)
    x = _cpp["rast.terrain"](x, vlist, neighbors, unit == "degrees", 0, opt)
    return messages(x, "terrain")


def sieve(
    x: SpatRaster,
    threshold: int,
    directions: int = 8,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Remove small patches — like R ``sieve()``."""
    opt = _opt(filename, **kw)
    x = _cpp["rast.sieve"](x, int(threshold), directions, opt)
    return messages(x, "sieve")


def rectify(
    x: SpatRaster,
    method: str = "bilinear",
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Rectify a rotated raster — like R ``rectify()``."""
    opt = _opt(filename, **kw)
    aoi = SpatRaster()
    x = x.rectify(method, aoi, 0, True, opt)
    return messages(x, "rectify")


def stretch(
    x: SpatRaster,
    minv: float = 0.0,
    maxv: float = 255.0,
    minq: float = 0.0,
    maxq: float = 1.0,
    smin: float = float("nan"),
    smax: float = float("nan"),
    bylayer: bool = True,
    maxcell: int = 500_000,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Stretch values — like R ``stretch()``."""
    opt = _opt(filename, **kw)
    x = _cpp["rast.stretch"](x, minv, maxv, minq, maxq, smin, smax, bylayer, maxcell, opt)
    return messages(x, "stretch")


def scale_linear(
    x: SpatRaster,
    min: float = 0.0,
    max: float = 1.0,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Linear rescale to [min, max] — like R ``scale_linear()``."""
    opt = _opt(filename, **kw)
    x = x.scale_linear(min, max, opt)
    return messages(x, "scale_linear")


def scale_raster(
    x: SpatRaster,
    center: Union[bool, List[float]] = True,
    scale: Union[bool, List[float]] = True,
) -> SpatRaster:
    """Center and/or scale — like R ``scale()``."""
    opt = SpatOptions()
    if isinstance(center, bool):
        docenter = center
        cv: List[float] = []
    else:
        docenter = True
        cv = [float(c) for c in center]
    if isinstance(scale, bool):
        doscale = scale
        sv: List[float] = []
    else:
        doscale = True
        sv = [float(s) for s in scale]
    x = x.scale(cv, docenter, sv, doscale, opt)
    return messages(x, "scale")


def quantile_raster(
    x: SpatRaster,
    probs: Optional[List[float]] = None,
    na_rm: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Per-cell quantiles — like R ``quantile()`` on SpatRaster."""
    if probs is None:
        probs = [0.0, 0.25, 0.5, 0.75, 1.0]
    opt = _opt(filename, **kw)
    x = x.quantile(probs, na_rm, opt)
    return messages(x, "quantile")


def atan_2(y: SpatRaster, x: SpatRaster, filename: str = "", **kw: Any) -> SpatRaster:
    """Arctan2 — like R ``atan_2()``."""
    opt = _opt(filename, **kw)
    y = y.atan2(x, opt)
    return messages(y, "atan_2")


# ── Raster processing ────────────────────────────────────────────────────────

def crop(
    x: SpatRaster,
    y: Any,
    snap: str = "near",
    mask: bool = False,
    touches: bool = True,
    extend: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """
    Crop a raster — like R ``crop()``.

    *y* can be a :class:`SpatExtent`, :class:`SpatRaster`, or :class:`SpatVector`.
    """
    opt = _opt(filename, **kw)
    if mask and isinstance(y, SpatVector):
        x = x.crop_mask(y, snap, touches, extend, opt)
        return messages(x, "crop")
    if isinstance(y, SpatRaster):
        e = y.extent
    elif isinstance(y, SpatVector):
        e = y.extent()
    elif isinstance(y, SpatExtent):
        e = y
    else:
        raise TypeError("crop: y must be SpatExtent, SpatRaster, or SpatVector")
    x = _cpp["rast.crop"](x, e, snap, extend, opt)
    return messages(x, "crop")


def mask(
    x: SpatRaster,
    mask_obj: Any,
    inverse: bool = False,
    mask_values: Optional[List[float]] = None,
    update_value: float = float("nan"),
    touches: bool = True,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """
    Mask a raster — like R ``mask()``.

    *mask_obj* can be a :class:`SpatRaster`, :class:`SpatVector`, or
    :class:`SpatExtent`.
    """
    opt = _opt(filename, **kw)
    if isinstance(mask_obj, SpatRaster):
        mv = mask_values if mask_values is not None else [float("nan")]
        x = x.mask_raster(mask_obj, inverse, mv, float(update_value), opt)
    elif isinstance(mask_obj, (SpatVector, SpatExtent)):
        if isinstance(mask_obj, SpatExtent):
            from .vect import vect as _vect
            from .crs import crs as _crs
            mask_obj = _vect(mask_obj, crs=_crs(x))
        x = x.mask_vector(mask_obj, inverse, float(update_value), touches, opt)
    else:
        raise TypeError("mask: mask_obj must be SpatRaster, SpatVector, or SpatExtent")
    return messages(x, "mask")


def project_raster(
    x: SpatRaster,
    y: Union[SpatRaster, str],
    method: str = "bilinear",
    mask: bool = False,
    align_only: bool = False,
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """
    Reproject a raster — like R ``project()`` on SpatRaster.

    *y* is either a target :class:`SpatRaster` (geometry template) or a CRS string.
    """
    opt = _opt(filename, **kw)
    if isinstance(y, SpatRaster):
        x = x.warp(y, "", method, mask, align_only, False, opt)
    else:
        x = x.warp(SpatRaster(), str(y), method, mask, align_only, False, opt)
    return messages(x, "project")


def resample(
    x: SpatRaster,
    y: SpatRaster,
    method: str = "bilinear",
    filename: str = "",
    **kw: Any,
) -> SpatRaster:
    """Resample to match *y* — like R ``resample()``."""
    opt = _opt(filename, **kw)
    x = x.warp(y, "", method, False, False, True, opt)
    return messages(x, "resample")


def intersect_rast(x: SpatRaster, y: SpatRaster) -> SpatRaster:
    """Intersection geometry and values — like R ``intersect()`` on SpatRaster."""
    opt = SpatOptions()
    x = x.intersect(y, opt)
    return messages(x, "intersect")


# ── Vector ───────────────────────────────────────────────────────────────────

def project_vector(
    x: SpatVector,
    y: Union[SpatRaster, SpatVector, str],
    partial: bool = False,
) -> SpatVector:
    """Reproject a vector — like R ``project()`` on SpatVector."""
    if not isinstance(y, str):
        if hasattr(y, "get_crs"):
            y = y.get_crs("wkt")
        else:
            y = str(y)
    x = x.project(y, partial)
    return messages(x, "project")


def shift_vect(x: SpatVector, dx: float = 0.0, dy: float = 0.0) -> SpatVector:
    """Shift vector coordinates — like R ``shift()`` on SpatVector."""
    x = x.shift(dx, dy)
    return messages(x, "shift")


def rotate_vect(x: SpatVector, longitude: float = 0.0, left: bool = True) -> SpatVector:
    """Rotate vector longitude — like R ``rotate()`` on SpatVector (non-split)."""
    x = x.rotate_longitude(longitude, left)
    return messages(x, "rotate")


def rescale_vect(
    x: SpatVector,
    fx: float = 0.5,
    fy: Optional[float] = None,
    x0: Optional[float] = None,
    y0: Optional[float] = None,
) -> SpatVector:
    """Rescale vector geometry — like R ``rescale()`` on SpatVector."""
    if fy is None:
        fy = fx
    e = x.extent().vector
    _x0 = x0 if x0 is not None else (e[0] + e[1]) / 2.0
    _y0 = y0 if y0 is not None else (e[2] + e[3]) / 2.0
    x = x.rescale(fx, fy, float(_x0), float(_y0))
    return messages(x, "rescale")


def trans_vect(x: SpatVector) -> SpatVector:
    """Transpose a vector — like R ``t()`` on SpatVector."""
    x = x.transpose()
    return messages(x, "trans")


# ── Scale-offset ─────────────────────────────────────────────────────────────

def scoff(x: SpatRaster) -> Any:
    """
    Return scale and offset per layer — like R ``scoff()``.

    Returns a list of ``[scale, offset]`` pairs.
    """
    so = x.getScaleOffset()
    scales = so[0]
    offsets = so[1]
    return [[s, o] for s, o in zip(scales, offsets)]


def scoff_set(
    x: SpatRaster,
    value: Optional[List[List[float]]] = None,
) -> SpatRaster:
    """
    Set scale and offset per layer — like R ``scoff(x) <- value``.

    *value* is a list of ``[scale, offset]`` pairs, one per layer.
    Use ``None`` to reset to ``[1, 0]`` for all layers.
    """
    x = x.deepcopy()
    n = x.nlyr()
    if value is None:
        x.setScaleOffset([1.0] * n, [0.0] * n)
    else:
        sc = [float(v[0]) if v[0] is not None else 1.0 for v in value]
        of = [float(v[1]) if v[1] is not None else 0.0 for v in value]
        x.setScaleOffset(sc, of)
        x.setValueType(0)
    return messages(x, "scoff<-")
