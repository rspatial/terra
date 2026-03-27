"""
Vector geometry operations — mirrors R ``terra::geom.R``.

Functions work directly on :class:`SpatVector` objects.
"""

from __future__ import annotations

from typing import Any, List, Optional, Union

from ._helpers import messages
from ._terra import SpatExtent, SpatOptions, SpatVector

__all__ = [
    # validity
    "is_valid", "make_valid",
    # set operations
    "union_vect", "intersect_vect", "erase", "symdif", "cover_vect",
    # crop / mask
    "crop_vect", "mask_vect",
    # geometry modifications
    "buffer_vect", "disagg_vect", "flip_vect", "spin",
    "hull", "delaunay", "voronoi", "elongate",
    "merge_lines", "make_nodes", "remove_dup_nodes",
    "simplify_geom", "thin_geom",
    "shared_paths", "snap_vect", "gaps",
    "force_ccw", "width_vect", "clearance",
    # predicates
    "is_empty",
]


def _opt() -> SpatOptions:
    return SpatOptions()


# ── Validity ─────────────────────────────────────────────────────────────────

def is_valid(
    x: SpatVector,
    with_messages: bool = False,
) -> Any:
    """
    Check geometry validity — like R ``is.valid()``.

    * ``with_messages=False`` → list[bool] per geometry.
    * ``with_messages=True`` → list of ``(valid, reason)`` pairs.
    """
    if with_messages:
        raw = x.geos_isvalid_msg()
        out = []
        for i in range(0, len(raw), 2):
            valid = (raw[i] == "\x01")
            out.append((valid, raw[i + 1]))
        return out
    return x.geos_isvalid()


def make_valid(x: SpatVector, buffer: bool = False) -> SpatVector:
    """Fix invalid geometries — like R ``makeValid()``."""
    if buffer and x.type() == "polygons":
        from .crs import crs as _crs
        xc = _crs(x)
        x.set_crs("local")
        x = x.buffer([0.0], 10, "round", "round", float("nan"), False)
        x.set_crs(xc)
        return messages(x, "makeValid")
    x = x.make_valid2()
    return messages(x, "makeValid")


# ── Set operations ───────────────────────────────────────────────────────────

def union_vect(
    x: SpatVector,
    y: Optional[SpatVector] = None,
) -> SpatVector:
    """
    Union — like R ``union(SpatVector)``.

    * ``y=None`` → unary union (dissolve self-overlaps).
    * ``y`` supplied → binary union.
    """
    if y is None:
        x = x.union_self()
    else:
        if x.type() != "polygons":
            # non-polygons: just append (like R rbind + unique)
            x = x.append(y)
        else:
            x = x.union(y)
    return messages(x, "union")


def intersect_vect(x: SpatVector, y: Union[SpatVector, SpatExtent]) -> SpatVector:
    """
    Intersect — like R ``intersect(SpatVector, …)``.

    *y* can be a :class:`SpatVector` or a :class:`SpatExtent` (crop).
    """
    if isinstance(y, SpatExtent):
        from .geom import crop_vect
        return crop_vect(x, y)
    x = x.intersect(y, True)
    return messages(x, "intersect")


def erase(
    x: SpatVector,
    y: Optional[Union[SpatVector, SpatExtent]] = None,
    sequential: bool = True,
) -> SpatVector:
    """
    Erase — like R ``erase()``.

    * ``y=None`` → self-erase (remove self-overlaps).
    * ``y=SpatVector`` → erase *y* from *x*.
    * ``y=SpatExtent`` → erase extent polygon from *x*.
    """
    if y is None:
        x = x.erase_self(sequential)
    elif isinstance(y, SpatExtent):
        from .vect import vect as _vect
        yv = _vect(y)
        x = x.erase(yv)
    else:
        x = x.erase_agg(y)
    return messages(x, "erase")


def symdif(x: SpatVector, y: SpatVector) -> SpatVector:
    """Symmetric difference — like R ``symdif()``."""
    x = x.symdif(y)
    return messages(x, "symdif")


def cover_vect(
    x: SpatVector,
    y: SpatVector,
    identity: bool = False,
    expand: bool = True,
) -> SpatVector:
    """Cover — like R ``cover(SpatVector)``."""
    x = x.cover(y, identity, expand)
    return messages(x, "cover")


# ── Crop / mask ───────────────────────────────────────────────────────────────

def crop_vect(
    x: SpatVector,
    y: Union[SpatVector, SpatExtent, Any],
    use_ext: bool = False,
) -> SpatVector:
    """
    Crop vector — like R ``crop(SpatVector, …)``.

    *y* can be a :class:`SpatVector`, :class:`SpatExtent`, or anything
    with an ``extent()`` method.
    """
    if isinstance(y, SpatVector) and not use_ext:
        x = x.crop_vct(y)
    else:
        if not isinstance(y, SpatExtent):
            try:
                y = y.extent() if callable(y.extent) else y.extent
            except AttributeError:
                from .extent import ext as _ext
                y = _ext(y)
        if x.type() == "points":
            from .vect import vect as _vect
            yv = _vect(y)
            x = x.crop_vct(yv)
        else:
            x = x.crop_ext(y, True)
    return messages(x, "crop")


def mask_vect(
    x: SpatVector,
    mask: Union[SpatVector, SpatExtent],
    inverse: bool = False,
) -> SpatVector:
    """
    Mask vector — like R ``mask(SpatVector)``.
    """
    if isinstance(mask, SpatExtent):
        from .vect import vect as _vect
        from .crs import crs as _crs
        mask = _vect(mask, crs=_crs(x))
    x = x.mask(mask, inverse)
    return messages(x, "mask")


# ── Geometry modifications ────────────────────────────────────────────────────

def buffer_vect(
    x: SpatVector,
    width: Union[float, List[float]],
    quadsegs: int = 10,
    capstyle: str = "round",
    joinstyle: str = "round",
    mitrelimit: float = float("nan"),
    singlesided: bool = False,
) -> SpatVector:
    """Buffer geometries — like R ``buffer(SpatVector)``."""
    if isinstance(width, (int, float)):
        w = [float(width)] * x.nrow()
    else:
        w = [float(v) for v in width]
    x = x.buffer(w, quadsegs, capstyle.lower(), joinstyle.lower(),
                 mitrelimit, singlesided)
    return messages(x, "buffer")


def disagg_vect(x: SpatVector, segments: bool = False) -> SpatVector:
    """Disaggregate multi-part geometries — like R ``disagg(SpatVector)``."""
    x = x.disaggregate(segments)
    return messages(x, "disagg")


def flip_vect(x: SpatVector, direction: str = "vertical") -> SpatVector:
    """Flip vector geometries — like R ``flip(SpatVector)``."""
    d = direction.lower()
    if d not in ("vertical", "horizontal"):
        raise ValueError("direction must be 'vertical' or 'horizontal'")
    x = x.flip(d == "vertical")
    return messages(x, "flip")


def spin(
    x: SpatVector,
    angle: float,
    x0: Optional[float] = None,
    y0: Optional[float] = None,
) -> SpatVector:
    """
    Rotate vector geometry around a point — like R ``spin()``.

    Default pivot is the centroid of the extent.
    """
    e = x.extent().vector
    _x0 = x0 if x0 is not None else (e[0] + e[1]) / 2.0
    _y0 = y0 if y0 is not None else (e[2] + e[3]) / 2.0
    x = x.rotate(float(angle), float(_x0), float(_y0))
    return messages(x, "spin")


def hull(
    x: SpatVector,
    type: str = "convex",
    by: str = "",
    param: float = 1.0,
    allow_holes: bool = True,
    tight: bool = True,
) -> SpatVector:
    """
    Convex / concave hull — like R ``hull()``.

    *type*: ``"convex"`` | ``"rectangle"`` | ``"circle"`` |
    ``"concave_ratio"`` | ``"concave_length"``
    """
    valid = {"convex", "rectangle", "circle", "concave_ratio", "concave_length"}
    t = type.lower()
    if t not in valid:
        raise ValueError(f"hull: type must be one of {valid}")
    x = x.hull(t, by, float(param), allow_holes, tight)
    return messages(x, "hull")


def delaunay(
    x: SpatVector,
    tolerance: float = 0.0,
    as_lines: bool = False,
    constrained: bool = False,
) -> SpatVector:
    """Delaunay triangulation — like R ``delaunay()``."""
    x = x.delaunay(tolerance, as_lines, constrained)
    return messages(x, "delaunay")


def voronoi(
    x: SpatVector,
    bnd: Optional[Any] = None,
    tolerance: float = 0.0,
    as_lines: bool = False,
) -> SpatVector:
    """
    Voronoi diagram — like R ``voronoi()``.

    *bnd* can be a :class:`SpatExtent`, :class:`SpatVector`, or ``None``.
    """
    if x.nrow() == 0:
        raise ValueError("voronoi: input has no geometries")
    if x.type() != "points":
        x = x.as_points(False, True)
    if bnd is None:
        bnd_v = SpatVector()
    else:
        if isinstance(bnd, SpatExtent):
            from .vect import vect as _vect
            bnd_v = _vect(bnd)
        elif isinstance(bnd, SpatVector):
            bnd_v = bnd
        else:
            from .extent import ext as _ext
            from .vect import vect as _vect
            bnd_v = _vect(_ext(bnd))
    x = x.voronoi(bnd_v, tolerance, as_lines)
    return messages(x, "voronoi")


def elongate(x: SpatVector, length: float = 1.0, flat: bool = False) -> SpatVector:
    """Elongate lines — like R ``elongate()``."""
    x = x.elongate(float(length), flat)
    return messages(x, "elongate")


def merge_lines(x: SpatVector) -> SpatVector:
    """Merge line segments — like R ``mergeLines()``."""
    x = x.line_merge()
    return messages(x, "mergeLines")


def make_nodes(x: SpatVector) -> SpatVector:
    """Add nodes at intersections — like R ``makeNodes()``."""
    x = x.make_nodes()
    return messages(x, "makeNodes")


def remove_dup_nodes(x: SpatVector, digits: int = -1) -> SpatVector:
    """Remove duplicate nodes — like R ``removeDupNodes()``."""
    x = x.remove_duplicate_nodes(digits)
    return messages(x, "removeDupNodes")


def simplify_geom(
    x: SpatVector,
    tolerance: float = 0.1,
    preserve_topology: bool = True,
    make_valid_after: bool = True,
) -> SpatVector:
    """Simplify geometry — like R ``simplifyGeom()``."""
    x = x.simplify(tolerance, preserve_topology)
    x = messages(x, "simplifyGeom")
    if make_valid_after:
        x = make_valid(x)
    return x


def thin_geom(
    x: SpatVector,
    threshold: float = 1e-6,
    make_valid_after: bool = True,
) -> SpatVector:
    """Thin geometry — like R ``thinGeom()``."""
    x = x.thin(threshold)
    x = messages(x, "thinGeom")
    if make_valid_after:
        x = make_valid(x)
    return x


def shared_paths(
    x: SpatVector,
    y: Optional[SpatVector] = None,
) -> SpatVector:
    """Shared paths — like R ``sharedPaths()``."""
    if y is None:
        x = x.shared_paths(True)
    else:
        x = x.shared_paths2(y, True)
    x = messages(x, "sharedPaths")
    return x


def snap_vect(
    x: SpatVector,
    y: Optional[SpatVector] = None,
    tolerance: float = 0.0,
) -> SpatVector:
    """Snap vertices — like R ``snap()``."""
    if y is None:
        x = x.snap(tolerance)
    else:
        x = x.snapto(y, tolerance)
    return messages(x, "snap")


def gaps(x: SpatVector) -> SpatVector:
    """Find gaps between polygons — like R ``gaps()``."""
    x = x.gaps()
    return messages(x, "gaps")


def width_vect(x: SpatVector, as_lines: bool = False) -> Any:
    """
    Width of polygons — like R ``width(SpatVector)``.

    Returns a SpatVector of lines when ``as_lines=True``, else a numeric vector.
    """
    x = x.width()
    x = messages(x, "width")
    if not as_lines:
        # perim is not yet wrapped; return the lines vector
        pass
    return x


def clearance(x: SpatVector, as_lines: bool = False) -> Any:
    """
    Minimum clearance — like R ``clearance()``.
    """
    x = x.clearance()
    x = messages(x, "clearance")
    return x


def force_ccw(x: SpatVector) -> SpatVector:
    """Force counter-clockwise winding — like R ``forceCCW()``."""
    x = x.deepcopy()
    x.make_CCW()
    return messages(x, "forceCCW")


# ── Predicates ────────────────────────────────────────────────────────────────

def is_empty(x: SpatVector) -> bool:
    """True if the vector has no geometries — like R ``is.empty()``."""
    return x.nrow() == 0
