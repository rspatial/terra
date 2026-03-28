"""Create SpatRasters — parallel to R ``terra::rast()``."""

from __future__ import annotations

from typing import Any, List, Optional, Sequence, Union

from ._helpers import character_crs, messages
from ._terra import SpatExtent, SpatOptions, SpatRaster, SpatVector

__all__ = ["rast"]


def _rast_from_file(
    path: str,
    *,
    subds: int = 0,
    drivers: Optional[List[str]] = None,
    opts: Optional[List[str]] = None,
    noflip: bool = False,
    guess_crs: bool = True,
    domains: Optional[List[str]] = None,
) -> SpatRaster:
    """Single-file raster using the C++ constructor (see R ``rast(character)``)."""
    f = [path]
    if drivers is None:
        drivers = []
    if opts is None:
        opts = []
    if domains is None:
        domains = []

    if subds < 1:
        subds_idx = -1
        subname: List[str] = []
    else:
        subds_idx = subds - 1
        subname = []

    r = SpatRaster(
        f,
        [subds_idx],
        subname,
        False,
        drivers,
        opts,
        [],
        noflip,
        guess_crs,
        domains,
    )
    return messages(r, "rast")


def rast(
    x: Any = None,
    *,
    nrows: int = 180,
    ncols: int = 360,
    nlyrs: int = 1,
    xmin: float = -180.0,
    xmax: float = 180.0,
    ymin: float = -90.0,
    ymax: float = 90.0,
    crs: Optional[str] = None,
    extent: Any = None,
    resolution: Optional[Union[float, Sequence[float]]] = None,
    vals: Any = None,
    names: Optional[Union[str, Sequence[str]]] = None,
    **kwargs: Any,
) -> SpatRaster:
    """
    Create a :class:`SpatRaster`, like R ``terra::rast()``.

    **From dimensions** (no positional *x*): same defaults as R ``rast()`` for
    a global lon/lat grid when CRS is omitted.

    **From a file path** (``x`` is ``str``): open with GDAL.

    **From a list of paths** (``x`` is a list of ``str``): combine sources.

    **From** :class:`SpatRaster`: return ``deepcopy``.

    **From** :class:`SpatExtent`: use extent as bounding box for a new empty raster.

    **From** :class:`SpatVector`: use ``x.extent()`` (and ``crs(x)`` if ``crs`` is
    omitted), like R ``rast(SpatVector, ...)``.

    Other keyword arguments are reserved for future parity with R ``rast(...)``.
    """
    del kwargs  # reserved

    # --- SpatRaster (copy) ---
    if isinstance(x, SpatRaster):
        return messages(x.deepcopy(), "rast")

    # --- From file path(s) ---
    if isinstance(x, str):
        return _rast_from_file(x)

    if isinstance(x, (list, tuple)) and len(x) > 0 and all(isinstance(s, str) for s in x):
        paths = list(x)
        if len(paths) == 1:
            return _rast_from_file(paths[0])
        out = _rast_from_file(paths[0])
        opt = SpatOptions()
        for i in range(1, len(paths)):
            r2 = _rast_from_file(paths[i])
            out.addSource(r2, True, opt)
        return messages(out, "rast")

    # --- SpatExtent ---
    if isinstance(x, SpatExtent):
        v = x.vector
        return rast(
            None,
            nrows=nrows,
            ncols=ncols,
            nlyrs=nlyrs,
            xmin=v[0],
            xmax=v[1],
            ymin=v[2],
            ymax=v[3],
            crs=crs,
            extent=None,
            resolution=resolution,
            vals=vals,
            names=names,
        )

    # --- SpatVector (R rast.R: extent + crs from vector) ---
    if isinstance(x, SpatVector):
        ext = x.extent()
        v = ext.vector
        vc = x.get_crs("wkt")
        crs_use = crs if crs is not None else (vc if vc else None)
        return rast(
            None,
            nrows=nrows,
            ncols=ncols,
            nlyrs=nlyrs,
            xmin=float(v[0]),
            xmax=float(v[1]),
            ymin=float(v[2]),
            ymax=float(v[3]),
            crs=crs_use,
            extent=None,
            resolution=resolution,
            vals=vals,
            names=names,
        )

    # --- new empty raster (x is None) ---
    if x is not None:
        raise TypeError(
            "rast: unsupported type for x; use None, str, list[str], SpatRaster, "
            "SpatExtent, or SpatVector"
        )

    nrows = int(round(nrows))
    ncols = int(round(ncols))
    nlyrs = int(round(nlyrs))
    if ncols < 1 or nrows < 1 or nlyrs < 1:
        raise ValueError("rast: nrows, ncols, nlyrs must be >= 1")

    if extent is not None:
        if hasattr(extent, "vector"):
            e = extent.vector
        else:
            e = list(extent)
        xmin, xmax, ymin, ymax = float(e[0]), float(e[1]), float(e[2]), float(e[3])
    else:
        xmin, xmax, ymin, ymax = float(xmin), float(xmax), float(ymin), float(ymax)

    if xmin >= xmax or ymin >= ymax:
        raise ValueError("rast: invalid extent")

    if crs is None:
        if xmin > -360.01 and xmax < 360.01 and ymin > -90.01 and ymax < 90.01:
            crs_use = "OGC:CRS84"
        else:
            crs_use = ""
    else:
        crs_use = character_crs(crs, "rast")

    r = SpatRaster(
        [max(1, nrows), max(1, ncols), max(1, nlyrs)],
        [xmin, xmax, ymin, ymax],
        crs_use,
    )
    r = messages(r, "rast")

    if resolution is not None:
        if isinstance(resolution, (int, float)):
            rx = ry = float(resolution)
        else:
            seq = list(resolution)
            if len(seq) < 2:
                raise ValueError("rast: resolution must be one or two numbers")
            rx, ry = float(seq[0]), float(seq[1])
        r = messages(r.set_resolution(rx, ry), "rast")

    if vals is not None:
        try:
            import numpy as np  # type: ignore

            if isinstance(vals, np.ndarray):
                flat = vals.astype("float64", copy=False).ravel(order="C")
                data = flat.tolist()
            else:
                data = list(vals)
        except ImportError:
            data = list(vals)

        opt = SpatOptions()
        r.setValues(data, opt)
        r = messages(r, "rast")

    if names is not None:
        if isinstance(names, str):
            nm = [names]
        else:
            nm = [str(s) for s in names]
        r.setNames(nm, True)

    return r
