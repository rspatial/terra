"""
spatvec.py — SpatVector geometry access and spatial measurements.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages, _getSpatDF


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# geometry type predicates
# ---------------------------------------------------------------------------

def geomtype(x: SpatVector) -> str:
    """Return the geometry type of *x*: ``'points'``, ``'lines'``, or ``'polygons'``."""
    # Must use C++ ``type()``, not ``geomtype()`` — the latter is patched onto
    # SpatVector by register_methods() and delegates back here (recursion).
    return x.type()


def is_lines(x: SpatVector) -> bool:
    """Return True if *x* is a line geometry."""
    return geomtype(x) == "lines"


def is_polygons(x: SpatVector) -> bool:
    """Return True if *x* is a polygon geometry."""
    return geomtype(x) == "polygons"


def is_points(x: SpatVector) -> bool:
    """Return True if *x* is a point geometry."""
    return "points" in geomtype(x)


# ---------------------------------------------------------------------------
# geometry access
# ---------------------------------------------------------------------------

def geom(
    x: SpatVector,
    wkt: bool = False,
    hex: bool = False,
    wkb: bool = False,
    as_df: bool = False,
    as_list: bool = False,
) -> Union[np.ndarray, List[str], List]:
    """
    Return the geometry of *x*.

    Parameters
    ----------
    x : SpatVector
    wkt : bool
        Return WKT strings.
    hex : bool
        Return hexadecimal WKB strings.
    wkb : bool
        Return raw WKB bytes.
    as_df : bool
        Return as a DataFrame (requires pandas).
    as_list : bool
        Return as a list of coordinate arrays.

    Returns
    -------
    numpy.ndarray, list of str, or list depending on arguments.
    """
    if hex:
        return list(x.hex())
    if wkt:
        return list(x.getGeometryWKT())
    if as_list:
        return x.get_geometryList("x", "y")
    if wkb:
        return x.wkb_raw()
    g = x.get_geometry()
    raw = list(g)
    n_cols = len(raw)
    if n_cols == 0:
        return np.empty((0, 5))
    max_len = max(len(col) for col in raw)
    arr = np.full((max_len, n_cols), float("nan"))
    for i, col in enumerate(raw):
        arr[:len(col), i] = col
    col_names = ["geom", "part", "x", "y", "hole"][:n_cols]
    if as_df:
        try:
            import pandas as pd
            return pd.DataFrame(arr, columns=col_names)
        except ImportError:
            pass
    return arr


def crds(
    x: Union[SpatVector, SpatRaster],
    as_df: bool = False,
    as_list: bool = False,
    na_rm: bool = True,
    na_all: bool = False,
) -> Union[np.ndarray, "pd.DataFrame", List]:
    """
    Return the coordinates of *x*.

    Parameters
    ----------
    x : SpatVector or SpatRaster
    as_df : bool
        Return as a DataFrame.
    as_list : bool
        Return as a list of arrays (SpatVector only).
    na_rm : bool
        Skip NA cells (SpatRaster only).
    na_all : bool
        Treat as NA only when all layers are NA (SpatRaster only).

    Returns
    -------
    numpy.ndarray, shape (n, 2), columns [x, y], or DataFrame/list.
    """
    if isinstance(x, SpatRaster):
        opt = _opt()
        out = x.crds(na_rm, na_all, opt)
        messages(x)
        raw = list(out)
        if len(raw) < 2:
            return np.empty((0, 2))
        arr = np.column_stack([raw[0], raw[1]])
    else:
        if as_list:
            gt = geomtype(x)
            if gt == "lines":
                return x.linesNA()
            elif gt == "polygons":
                return x.polygonsList()
            else:
                return x.coordinates()
        raw = x.coordinates()
        raw_list = list(raw)
        if len(raw_list) < 2:
            return np.empty((0, 2))
        arr = np.column_stack([raw_list[0], raw_list[1]])

    if as_df:
        try:
            import pandas as pd
            return pd.DataFrame(arr, columns=["x", "y"])
        except ImportError:
            pass
    return arr


# ---------------------------------------------------------------------------
# spatial measurements
# ---------------------------------------------------------------------------

def expanse(
    x: SpatVector,
    unit: str = "m",
    transform: bool = True,
) -> np.ndarray:
    """
    Compute the area of each feature.

    Parameters
    ----------
    x : SpatVector
        Must contain polygon geometries.
    unit : str
        ``"m"`` (square metres, default) or ``"km"`` (square kilometres).
    transform : bool
        Transform geographic coordinates to an equal-area projection.

    Returns
    -------
    numpy.ndarray of float (one value per feature).
    """
    a = x.area(unit, transform, [])
    messages(x, "expanse")
    return np.abs(np.array(a, dtype=float))


def perim(x: SpatVector) -> np.ndarray:
    """
    Compute the perimeter (or length) of each feature.

    Parameters
    ----------
    x : SpatVector

    Returns
    -------
    numpy.ndarray of float.
    """
    p = x.length()
    messages(x, "perim")
    return np.array(p, dtype=float)


def nseg(x: SpatVector) -> np.ndarray:
    """Return the number of segments per feature."""
    s = x.nsegments()
    messages(x, "nseg")
    return np.array(s, dtype=int)


# ---------------------------------------------------------------------------
# fill holes
# ---------------------------------------------------------------------------

def fill_holes(x: SpatVector, inverse: bool = False) -> SpatVector:
    """
    Fill interior holes in polygon geometries.

    Parameters
    ----------
    x : SpatVector
    inverse : bool
        If True, return only the holes (i.e. extract holes).

    Returns
    -------
    SpatVector
    """
    if inverse:
        xc = x.get_holes()
    else:
        xc = x.remove_holes()
    return messages(xc, "fillHoles")


# ---------------------------------------------------------------------------
# as_data_frame for SpatVector
# ---------------------------------------------------------------------------

def vect_as_df(
    x: SpatVector,
    geom: Optional[str] = None,
) -> "pd.DataFrame":
    """
    Convert the attribute table of *x* to a pandas DataFrame.

    Parameters
    ----------
    x : SpatVector
    geom : str, optional
        If provided, append geometry: ``"WKT"``, ``"HEX"``, or ``"XY"``
        (points only).

    Returns
    -------
    pandas.DataFrame
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for vect_as_df()")

    from ._helpers import _getSpatDF
    d = _getSpatDF(x.df)
    if d is None:
        d = pd.DataFrame()

    if geom is not None:
        geom = geom.upper()
        if geom == "XY":
            if not is_points(x):
                raise ValueError('geom="XY" is only valid for point geometries')
            xy = crds(x, as_df=True)
            d = pd.concat([d.reset_index(drop=True), xy.reset_index(drop=True)], axis=1)
        elif geom in ("WKT", "HEX"):
            d = d.copy()
            d["geometry"] = geom_as_wkt(x) if geom == "WKT" else list(x.hex())
    return d


def geom_as_wkt(x: SpatVector) -> List[str]:
    """Return WKT representation for each feature."""
    return list(x.getGeometryWKT())
