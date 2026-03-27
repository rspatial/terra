"""
Vector geometry operations.

Functions operate directly on :class:`SpatVector` objects and cover set
operations (union, intersect, erase, cover, symmetric difference), geometry
modification (buffer, hull, simplify, snap, flip, rotate, …), topology
utilities, and predicates.
"""

from __future__ import annotations

from typing import Any, List, Optional, Union

from ._helpers import messages
from ._terra import SpatExtent, SpatOptions, SpatVector

# Captured before monkey-patching in methods.py
_cpp_vect_union      = SpatVector.union
_cpp_vect_intersect  = SpatVector.intersect
_cpp_vect_erase      = SpatVector.erase

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
    Check the validity of polygon geometries.

    Args:
        x: SpatVector of polygons (or any geometry type; non-polygon
            geometries are always considered valid).
        with_messages: If True, return a list of ``(valid, reason)`` tuples
            with one entry per geometry, where ``reason`` is a string
            describing the invalidity (empty string for valid geometries).
            If False (default), return a list of booleans.

    Returns:
        List of bool, or list of (bool, str) tuples when
        ``with_messages=True``.
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
    """
    Attempt to fix invalid polygon geometries.

    Args:
        x: SpatVector of polygons.
        buffer: If True, a zero-width buffer is used to create valid polygons.
            Use with caution: this method may result in data loss, for example
            only a single part of a self-intersecting polygon may be
            preserved.

    Returns:
        SpatVector with valid geometries.
    """
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
    Compute the union of vector geometries.

    When ``y`` is provided, the two SpatVectors are combined and overlapping
    polygons between (not within) objects are intersected. For lines and
    points the datasets are simply combined (equivalent to ``rbind``).
    Attributes are joined.

    When ``y`` is omitted, overlapping polygons within ``x`` are
    intersected (self-union). Original attributes are lost; new attributes
    record how many and which input polygons overlapped.

    Args:
        x: SpatVector.
        y: SpatVector to union with ``x``. If not provided, a self-union is
            performed.

    Returns:
        SpatVector.
    """
    if y is None:
        x = x.union_self()
    else:
        if x.type() != "polygons":
            x = x.append(y)
        else:
            x = _cpp_vect_union(x, y)
    return messages(x, "union")


def intersect_vect(x: SpatVector, y: Union[SpatVector, SpatExtent]) -> SpatVector:
    """
    Intersect a SpatVector with another SpatVector or a SpatExtent.

    When ``y`` is a SpatExtent, the operation is equivalent to
    :func:`crop_vect` with the extent.

    Intersecting points with points uses the extent of ``y``.  Intersecting
    points and lines is not supported; create a polygon buffer from the lines
    first and use that with ``intersect_vect``.

    Args:
        x: SpatVector.
        y: SpatVector or SpatExtent to intersect with.

    Returns:
        SpatVector containing only the overlapping parts of ``x`` and ``y``,
        with combined attributes.
    """
    if isinstance(y, SpatExtent):
        return crop_vect(x, y)
    x = _cpp_vect_intersect(x, y, True)
    return messages(x, "intersect")


def erase(
    x: SpatVector,
    y: Optional[Union[SpatVector, SpatExtent]] = None,
    sequential: bool = True,
) -> SpatVector:
    """
    Erase parts of a SpatVector.

    Three forms are supported:

    * **Self-erase** (``y=None``): erase the overlapping areas within ``x``
      itself. If ``sequential=True`` (default), each polygon is erased
      against the preceding ones in order. If ``sequential=False``, all
      overlapping areas are removed and only the areas covered by a single
      geometry are kept.
    * **Erase with SpatVector**: erase the areas of ``y`` from ``x``.
    * **Erase with SpatExtent**: erase the rectangular region from ``x``.

    Args:
        x: SpatVector.
        y: SpatVector, SpatExtent, or None (self-erase).
        sequential: Only used for self-erase. If True, polygons are erased
            sequentially. If False, all overlapping areas are erased at once.

    Returns:
        SpatVector.
    """
    if y is None:
        x = x.erase_self(sequential)
    elif isinstance(y, SpatExtent):
        from .vect import vect as _vect
        yv = _vect(y)
        x = _cpp_vect_erase(x, yv)
    else:
        x = x.erase_agg(y)
    return messages(x, "erase")


def symdif(x: SpatVector, y: SpatVector) -> SpatVector:
    """
    Compute the symmetrical difference of two polygon SpatVectors.

    The result contains the areas that are in ``x`` or ``y`` but not in both.

    Args:
        x: SpatVector of polygons.
        y: SpatVector of polygons.

    Returns:
        SpatVector.
    """
    x = x.symdif(y)
    return messages(x, "symdif")


def cover_vect(
    x: SpatVector,
    y: SpatVector,
    identity: bool = False,
    expand: bool = True,
) -> SpatVector:
    """
    Replace parts of ``x`` that overlap with ``y``.

    Areas of ``x`` that overlap with ``y`` are replaced by ``y``, or, if
    ``identity=True``, intersected with ``y``.

    Args:
        x: SpatVector.
        y: SpatVector used to cover or intersect overlapping areas of ``x``.
        identity: If True, overlapping areas are intersected rather than
            simply replaced.
        expand: If True (default), parts of ``y`` that are outside ``x`` are
            included in the result.

    Returns:
        SpatVector.
    """
    x = x.cover(y, identity, expand)
    return messages(x, "cover")


# ── Crop / mask ───────────────────────────────────────────────────────────────

def crop_vect(
    x: SpatVector,
    y: Union[SpatVector, SpatExtent, Any],
    use_ext: bool = False,
) -> SpatVector:
    """
    Cut out a geographic subset of a SpatVector.

    When cropping a SpatVector with another SpatVector, the minimum convex
    hull of ``y`` is used unless ``use_ext=True``. Unlike :func:`intersect_vect`,
    the attributes of ``y`` are not transferred to the result.

    You can crop a SpatVector with a rectangle (SpatExtent or any object with
    an extent).

    Args:
        x: SpatVector.
        y: SpatVector, SpatExtent, or any object that has an extent.
        use_ext: If True, use the bounding extent of ``y`` rather than ``y``
            itself as the crop region.

    Returns:
        SpatVector cropped to the region defined by ``y``.
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
    Select geometries of ``x`` that intersect (or do not intersect) with ``mask``.

    Args:
        x: SpatVector.
        mask: SpatVector or SpatExtent defining the mask region.
        inverse: If True, select geometries that do *not* intersect with
            ``mask``.

    Returns:
        SpatVector with a subset of the geometries of ``x``.
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
    """
    Calculate a buffer around vector geometries.

    The unit of ``width`` is metres when the CRS is longitude/latitude, and
    in the linear unit of the coordinate reference system otherwise (typically
    also metres).

    Args:
        x: SpatVector.
        width: Buffer width. A single value is applied to all geometries; a
            list of values applies a different width to each geometry.
        quadsegs: Number of line segments used to approximate a quarter circle
            (higher values produce smoother circles). Positive integer.
        capstyle: End cap style for lines. One of ``"round"``, ``"square"``,
            or ``"flat"``. Ignored for longitude/latitude CRS.
        joinstyle: Join style. One of ``"round"``, ``"mitre"``, or
            ``"bevel"``. Ignored for longitude/latitude CRS.
        mitrelimit: Upper bound on mitre joins, to prevent them from
            extending very far at acute angles. Ignored for
            longitude/latitude CRS.
        singlesided: If True, construct a buffer on only one side of each
            input line. Ignored for longitude/latitude CRS.

    Returns:
        SpatVector with buffered geometries.
    """
    if isinstance(width, (int, float)):
        w = [float(width)] * x.nrow()
    else:
        w = [float(v) for v in width]
    x = x.buffer(w, quadsegs, capstyle.lower(), joinstyle.lower(),
                 mitrelimit, singlesided)
    return messages(x, "buffer")


def disagg_vect(x: SpatVector, segments: bool = False) -> SpatVector:
    """
    Separate multi-part geometries into single-part geometries.

    For lines or polygons, setting ``segments=True`` further breaks each
    geometry into its individual line segments.

    Args:
        x: SpatVector.
        segments: If True, (poly-)lines and polygons are disaggregated into
            individual line segments.

    Returns:
        SpatVector with single-part geometries.
    """
    x = x.disaggregate(segments)
    return messages(x, "disagg")


def flip_vect(x: SpatVector, direction: str = "vertical") -> SpatVector:
    """
    Flip vector geometries along a vertical or horizontal axis.

    Args:
        x: SpatVector.
        direction: ``"vertical"`` to flip by rows (invert the y-axis) or
            ``"horizontal"`` to flip by columns (invert the x-axis).

    Returns:
        SpatVector.
    """
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
    Rotate a SpatVector around a pivot point.

    Args:
        x: SpatVector.
        angle: Rotation angle in degrees (counter-clockwise positive).
        x0: X coordinate of the pivot point. Defaults to the centre of the
            bounding extent.
        y0: Y coordinate of the pivot point. Defaults to the centre of the
            bounding extent.

    Returns:
        SpatVector.
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
    Compute a hull around SpatVector geometries.

    Args:
        x: SpatVector.
        type: Hull type. One of:

            * ``"convex"``         — convex hull (default)
            * ``"rectangle"``      — minimum bounding rotated rectangle
            * ``"circle"``         — minimum bounding circle
            * ``"concave_ratio"``  — concave hull specified by edge-length ratio
            * ``"concave_length"`` — concave hull specified by maximum edge length

        by: Name of an attribute variable in ``x`` to compute a separate hull
            for each group of geometries.
        param: Concaveness parameter for ``"concave_*"`` types. For
            ``"concave_ratio"``, a value between 0 and 1 (1 = convex hull,
            0 = maximum concaveness). For ``"concave_length"``, the maximum
            edge length (larger = more convex).
        allow_holes: If True, the output polygons may contain holes. Only
            used for ``"concave_*"`` types.
        tight: If True, the hull closely follows the outer boundaries of the
            input polygons. Only used for ``"concave_length"`` with polygon
            input.

    Returns:
        SpatVector of polygons.
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
    """
    Compute Delaunay triangles for point, line, or polygon nodes.

    Args:
        x: SpatVector.
        tolerance: Snapping tolerance (≥ 0). A value of 0 disables snapping.
        as_lines: If True, return triangles as lines without the outer
            boundary.
        constrained: If True, compute a constrained Delaunay triangulation.

    Returns:
        SpatVector of polygons (or lines when ``as_lines=True``).
    """
    x = x.delaunay(tolerance, as_lines, constrained)
    return messages(x, "delaunay")


def voronoi(
    x: SpatVector,
    bnd: Optional[Any] = None,
    tolerance: float = 0.0,
    as_lines: bool = False,
) -> SpatVector:
    """
    Compute the Voronoi diagram for points or the nodes of geometries.

    Args:
        x: SpatVector. Only point coordinates are used; if ``x`` contains
            lines or polygons their nodes are used.
        bnd: SpatVector or SpatExtent defining the outer boundary of the
            diagram. If None (default), the boundary is derived from ``x``.
        tolerance: Snapping tolerance (≥ 0). A value of 0 disables snapping.
        as_lines: If True, return Voronoi edges as lines without the outer
            boundary polygon.

    Returns:
        SpatVector of polygons (or lines when ``as_lines=True``).
    """
    if x.nrow() == 0:
        raise ValueError("voronoi: input has no geometries")
    if x.type() != "points":
        x = x.as_points(False, True)
    if bnd is None:
        bnd_v = SpatVector()
    elif isinstance(bnd, SpatExtent):
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
    """
    Extend lines by a given length at each end.

    Args:
        x: SpatVector of lines.
        length: Amount to elongate each end of each line. The unit is metres
            when the CRS is longitude/latitude; otherwise the same as the
            linear unit of the CRS (typically also metres).
        flat: If True, the earth's curvature is ignored for longitude/latitude
            data, and the unit is degrees rather than metres.

    Returns:
        SpatVector with elongated lines.
    """
    x = x.elongate(float(length), flat)
    return messages(x, "elongate")


def merge_lines(x: SpatVector) -> SpatVector:
    """
    Merge line segments that share endpoints into longer lines.

    Can also be used to convert closed line rings into polygons.

    Args:
        x: SpatVector of lines.

    Returns:
        SpatVector with merged lines.
    """
    x = x.line_merge()
    return messages(x, "mergeLines")


def make_nodes(x: SpatVector) -> SpatVector:
    """
    Create nodes (split lines) at every intersection.

    Args:
        x: SpatVector of lines or polygons.

    Returns:
        SpatVector with added nodes at all line intersections.
    """
    x = x.make_nodes()
    return messages(x, "makeNodes")


def remove_dup_nodes(x: SpatVector, digits: int = -1) -> SpatVector:
    """
    Remove duplicate nodes and optionally round coordinates.

    Args:
        x: SpatVector of lines or polygons.
        digits: Number of decimal digits used in rounding. If negative
            (default), no rounding is applied.

    Returns:
        SpatVector with duplicate nodes removed.
    """
    x = x.remove_duplicate_nodes(digits)
    return messages(x, "removeDupNodes")


def simplify_geom(
    x: SpatVector,
    tolerance: float = 0.1,
    preserve_topology: bool = True,
    make_valid_after: bool = True,
) -> SpatVector:
    """
    Reduce the number of nodes used to represent geometries.

    Args:
        x: SpatVector of lines or polygons.
        tolerance: Minimum distance between remaining nodes, in the units of
            the CRS (degrees for longitude/latitude).
        preserve_topology: If True, the topological relationships of the
            output geometries are preserved.
        make_valid_after: If True, :func:`make_valid` is run on the result
            to ensure that all output polygons are valid.

    Returns:
        SpatVector with simplified geometries.
    """
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
    """
    Thin geometries by removing nodes that are closer than a threshold.

    Args:
        x: SpatVector of lines or polygons.
        threshold: Minimum distance between nodes. Nodes closer than this
            value are removed.
        make_valid_after: If True, :func:`make_valid` is run on the result
            to ensure validity.

    Returns:
        SpatVector with thinned geometries.
    """
    x = x.thin(threshold)
    x = messages(x, "thinGeom")
    if make_valid_after:
        x = make_valid(x)
    return x


def shared_paths(
    x: SpatVector,
    y: Optional[SpatVector] = None,
) -> SpatVector:
    """
    Find the shared paths between line or polygon geometries.

    Args:
        x: SpatVector of lines or polygons.
        y: Second SpatVector to compare against. If None (default), shared
            paths are found among the geometries within ``x`` itself.

    Returns:
        SpatVector of lines with attributes ``id1`` and ``id2`` identifying
        which input geometries share each path.
    """
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
    """
    Snap the vertices of ``x`` to those of ``y`` (or to each other) within
    a distance tolerance.

    This makes nearby boundaries of geometries identical, resolving
    near-coincident edges.

    Args:
        x: SpatVector of lines or polygons.
        y: Target SpatVector to snap towards. If None (default), ``x`` is
            snapped to its own geometries.
        tolerance: Snapping tolerance — the maximum distance between vertices
            that will be snapped together.

    Returns:
        SpatVector with snapped vertices.
    """
    if y is None:
        x = x.snap(tolerance)
    else:
        x = x.snapto(y, tolerance)
    return messages(x, "snap")


def gaps(x: SpatVector) -> SpatVector:
    """
    Find the gaps (holes) between polygons.

    Args:
        x: SpatVector of polygons.

    Returns:
        SpatVector of polygons representing the gaps between the input
        polygons.
    """
    x = x.gaps()
    return messages(x, "gaps")


def width_vect(x: SpatVector, as_lines: bool = False) -> Any:
    """
    Compute the minimum diameter of each geometry.

    The minimum diameter is the width of the smallest band (strip defined by
    two parallel lines) that fully contains the geometry. This can be thought
    of as the smallest hole that the geometry can pass through with a single
    rotation.

    Args:
        x: SpatVector of lines or polygons.
        as_lines: If True, return the lines that define the minimum diameter
            as a SpatVector. If False (default), return numeric widths.

    Returns:
        Numeric values (one per geometry), or SpatVector of lines when
        ``as_lines=True``.
    """
    x = x.width()
    x = messages(x, "width")
    return x


def clearance(x: SpatVector, as_lines: bool = False) -> Any:
    """
    Compute the minimum clearance of each geometry.

    The minimum clearance is the smallest amount by which a vertex could be
    moved to produce an invalid polygon, a non-simple linestring, or a
    multipoint with repeated points. If the minimum clearance cannot be
    defined (e.g. for a single point), NA is returned.

    Args:
        x: SpatVector of lines or polygons.
        as_lines: If True, return the lines that define the clearance as a
            SpatVector. If False (default), return numeric clearance values.

    Returns:
        Numeric values (one per geometry), or SpatVector of lines when
        ``as_lines=True``.
    """
    x = x.clearance()
    x = messages(x, "clearance")
    return x


def force_ccw(x: SpatVector) -> SpatVector:
    """
    Force the outer rings of polygons to follow counter-clockwise winding order.

    Args:
        x: SpatVector of polygons.

    Returns:
        SpatVector with counter-clockwise polygon rings.
    """
    x = x.deepcopy()
    x.make_CCW()
    return messages(x, "forceCCW")


# ── Predicates ────────────────────────────────────────────────────────────────

def is_empty(x: SpatVector) -> bool:
    """
    Test whether a SpatVector contains no geometries.

    Args:
        x: SpatVector.

    Returns:
        True if ``x`` has no geometries, False otherwise.
    """
    return x.nrow() == 0
