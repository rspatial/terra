"""
distance.py — distance and buffer operations for SpatRaster and SpatVector.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages, spatoptions


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# raster buffer
# ---------------------------------------------------------------------------

def buffer_rast(
    x: SpatRaster,
    width: float,
    background: float = 0.0,
    include: bool = True,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Create a buffer around non-NA or non-zero cells.

    Parameters
    ----------
    x : SpatRaster
    width : float
        Buffer width in map units (or meters for lonlat).
    background : float
        Value for non-buffered cells.
    include : bool
        Include the source cells in the buffer.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    opt = spatoptions(filename, overwrite)
    xc = x.buffer(width, background, include, opt)
    return messages(xc, "buffer")


# ---------------------------------------------------------------------------
# raster distance
# ---------------------------------------------------------------------------

def distance_rast(
    x: SpatRaster,
    y: Optional[SpatVector] = None,
    *,
    target: float = float("nan"),
    exclude: Optional[float] = None,
    unit: str = "m",
    method: str = "haversine",
    maxdist: float = float("nan"),
    values: bool = False,
    rasterize: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Compute the distance from each cell to the nearest non-NA (or target) cell.

    Parameters
    ----------
    x : SpatRaster
    y : SpatVector, optional
        If provided, compute distance from *x* cells to *y* features.
    target : float
        Cell value to use as source (default: non-NA).
    exclude : float, optional
        Cell value to exclude from sources.
    unit : str
        ``"m"`` (metres) or ``"km"`` (kilometres).
    method : str
        ``"haversine"`` (default), ``"cosine"``, or ``"geo"``.
    maxdist : float
        Maximum distance to compute; cells beyond this get NA.
    values : bool
        Return the value at the nearest source instead of the distance.
    rasterize : bool
        When *y* is provided, rasterize *y* first (faster for many features).
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    method = method.lower()
    if method not in ("cosine", "haversine", "geo"):
        raise ValueError(f"method must be 'cosine', 'haversine', or 'geo'; got {method!r}")
    opt = spatoptions(filename, overwrite)
    if y is not None:
        xc = x.vectDistance(y, rasterize, unit, method, opt)
    else:
        keep_na = False
        exclude_val = float("nan")
        if exclude is not None:
            exclude_val = float(exclude)
            if np.isnan(exclude_val) and np.isnan(target):
                raise ValueError("'target' and 'exclude' must be different")
            if np.isnan(exclude_val):
                keep_na = True
        xc = x.rastDistance(target, exclude_val, keep_na, unit, True, method, values, maxdist, opt)
    return messages(xc, "distance")


def cost_dist(
    x: SpatRaster,
    target: float = 0.0,
    scale: float = 1.0,
    maxiter: int = 50,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Least-cost distance based on cell values as conductance weights.

    Parameters
    ----------
    x : SpatRaster
        Conductance surface (NA = barrier).
    target : float
        Value identifying source cells.
    scale : float
        Scale factor applied to cell values.
    maxiter : int
        Maximum iterations.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    opt = spatoptions(filename, overwrite)
    xc = x.costDistance(target, scale, max(maxiter, 2), False, opt)
    return messages(xc, "costDist")


def grid_dist(
    x: SpatRaster,
    target: Optional[float] = 0.0,
    scale: float = 1.0,
    maxiter: int = 50,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Grid (cumulative) distance from source cells.

    Parameters
    ----------
    x : SpatRaster
    target : float or None
        Source cell value (None → NA distance field).
    scale : float
    maxiter : int
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    opt = spatoptions(filename, overwrite)
    if target is None or np.isnan(float(target)):
        xc = x.gridDistance(float(scale), opt)
    else:
        xc = x.costDistance(float(target), float(scale), max(maxiter, 2), True, opt)
    return messages(xc, "gridDist")


# ---------------------------------------------------------------------------
# vector distance
# ---------------------------------------------------------------------------

def _test_for_lonlat(xy: np.ndarray) -> bool:
    """R ``test.for.lonlat`` — ``R/distance.R`` (used when *lonlat* is omitted)."""
    xr = (float(np.nanmin(xy[:, 0])), float(np.nanmax(xy[:, 0])))
    yr = (float(np.nanmin(xy[:, 1])), float(np.nanmax(xy[:, 1])))
    return bool(
        xr[0] >= -180.0
        and xr[1] <= 180.0
        and yr[0] > -90.0
        and yr[1] < 90.0
    )


def distance_xy(
    x: np.ndarray,
    *,
    lonlat: Optional[bool] = None,
    sequential: bool = False,
    pairs: bool = False,
    symmetrical: bool = True,
    unit: str = "m",
    method: str = "geo",
    use_nodes: bool = False,
) -> np.ndarray:
    """
    Mirror R ``distance(matrix, y=missing, ...)`` — ``R/distance.R`` lines 221–231.

    Builds CRS with ``+proj=utm +zone=1 +datum=WGS84`` when ``lonlat`` is
    false, else ``+proj=longlat +datum=WGS84``, then ``vect(matrix, crs=…)``
    and ``distance()`` on the vector (same defaults as R: ``method="geo"``).
    """
    from .vect import vect

    m = np.asarray(x, dtype=float)
    if m.ndim != 2 or m.shape[1] != 2:
        raise ValueError("distance_xy: x must be an array of shape (n, 2)")
    if lonlat is None:
        lonlat = _test_for_lonlat(m)
    crs = (
        "+proj=longlat +datum=WGS84"
        if lonlat
        else "+proj=utm +zone=1 +datum=WGS84"
    )
    v = vect(m, crs=crs)
    return distance_vect_self(
        v,
        sequential=sequential,
        pairs=pairs,
        symmetrical=symmetrical,
        unit=unit,
        method=method,
        use_nodes=use_nodes,
    )


def distance_vect_self(
    x: SpatVector,
    sequential: bool = False,
    pairs: bool = False,
    symmetrical: bool = True,
    unit: str = "m",
    method: str = "haversine",
    use_nodes: bool = False,
) -> np.ndarray:
    """
    Compute pairwise distances between all features of *x*.

    The underlying C++ implementation requires a defined CRS on *x*; otherwise
    it returns an empty result.

    **Method names** match R ``distance(SpatVector, …)``: only ``"geo"``,
    ``"haversine"``, and ``"cosine"`` are valid (R ``match.arg``).  There is
    no ``"euclidean"`` label: for a **projected** CRS, planar distances use the
    same C++ path regardless of which of the three names you pass.

    Parameters
    ----------
    x : SpatVector
    sequential : bool
        If True, return only between consecutive features.
    pairs : bool
        If True, return a 3-column matrix [from, to, distance].
    symmetrical : bool
        If True (and pairs), return only the upper triangle.
    unit : str
    method : str
    use_nodes : bool
        Use node coordinates for lines/polygons.

    Returns
    -------
    numpy.ndarray  (condensed distance vector, or matrix, or pair matrix).
    """
    method = method.lower()
    if method not in ("cosine", "haversine", "geo"):
        raise ValueError(
            "distance_vect_self: method must be 'geo', 'haversine', or 'cosine' "
            f"(same as R distance(SpatVector)); got {method!r}"
        )
    opt = _opt()
    d = x.distance_self(sequential, unit, method, use_nodes, opt)
    d_arr = np.array(d, dtype=float)
    if sequential:
        return d_arr
    n = x.nrow()
    need = n * (n - 1) // 2
    mat = np.full((n, n), float("nan"))
    idx = np.triu_indices(n, k=1)
    if d_arr.size != need:
        raise ValueError(
            f"distance_self: expected {need} pairwise distances, got {d_arr.size}"
        )
    mat[idx] = d_arr
    mat = mat + mat.T
    np.fill_diagonal(mat, 0.0)
    if pairs:
        if symmetrical:
            r, c = idx
        else:
            r, c = np.indices((n, n))
            r, c = r.ravel(), c.ravel()
            mask = r != c
            r, c = r[mask], c[mask]
        out = np.column_stack([r + 1, c + 1, mat[r, c]])
        return out
    return mat


def distance_vect(
    x: SpatVector,
    y: SpatVector,
    pairwise: bool = False,
    unit: str = "m",
    method: str = "haversine",
    use_nodes: bool = False,
) -> np.ndarray:
    """
    Compute distances from each feature in *x* to each feature in *y*.

    Parameters
    ----------
    x : SpatVector
    y : SpatVector
    pairwise : bool
        If True, compute n × 1 pairwise distances (requires nrow(x) == nrow(y)).
    unit : str
    method : str
    use_nodes : bool

    Returns
    -------
    numpy.ndarray, shape (nrow(x), nrow(y)) or (nrow(x),) if pairwise.
    """
    method = method.lower()
    opt = _opt()
    d = x.distance_other(y, pairwise, unit, method, use_nodes, opt)
    d_arr = np.array(d, dtype=float)
    if pairwise:
        return d_arr
    return d_arr.reshape(x.nrow(), y.nrow())


def distance_points(
    x: np.ndarray,
    y: np.ndarray,
    lonlat: Optional[bool] = None,
    pairwise: bool = False,
    unit: str = "m",
    method: str = "geo",
) -> np.ndarray:
    """
    Compute distances between two sets of coordinate pairs.

    Parameters
    ----------
    x : array-like, shape (n, 2)  [lon/x, lat/y]
    y : array-like, shape (m, 2)
    lonlat : bool, optional
        Whether coordinates are geographic.  Auto-detected if None.
    pairwise : bool
        If True, compute n pairwise distances (requires n == m).
    unit : str
        ``"m"`` or ``"km"``.
    method : str

    Returns
    -------
    numpy.ndarray, shape (n, m) or (n,) if pairwise.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if lonlat is None:
        def _looks_lonlat(arr):
            return (arr[:, 0].min() >= -180 and arr[:, 0].max() <= 180
                    and arr[:, 1].min() > -90 and arr[:, 1].max() < 90)
        lonlat = _looks_lonlat(x_arr) and _looks_lonlat(y_arr)

    m_scale = 1.0 if unit == "m" else 0.001
    v = SpatVector()
    d = v.point_distance(
        x_arr[:, 0].tolist(), x_arr[:, 1].tolist(),
        y_arr[:, 0].tolist(), y_arr[:, 1].tolist(),
        pairwise, m_scale, lonlat, method,
    )
    d_arr = np.array(d, dtype=float)
    if pairwise:
        return d_arr
    return d_arr.reshape(len(x_arr), len(y_arr))
