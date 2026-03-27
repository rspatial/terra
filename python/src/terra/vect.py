"""Create vectors — parallel to R ``terra::vect()``."""

from __future__ import annotations

import os
from typing import Any, List, Optional

from ._helpers import character_crs, messages
from ._terra import SpatExtent, SpatVector

__all__ = ["vect"]


def _looks_like_wkt(s: str) -> bool:
    # Match R's first-five-character test (POINT, MULTI, LINES, POLYG, EMPTY)
    t = s.strip().upper()[:5]
    return t in ("POINT", "MULTI", "LINES", "POLYG", "EMPTY") or s.strip().startswith("{")


def _normalize_path(path: str) -> str:
    p = path.strip()
    if p.startswith("http") and (p.endswith(".shp") or p.endswith(".gpkg")):
        return "/vsicurl/" + p
    if p.startswith("s3://"):
        return "/vsis3/" + p[5:]
    if os.path.isfile(p):
        try:
            return os.path.abspath(p)
        except OSError:
            return p
    return p


def vect(
    x: Any = None,
    *,
    layer: str = "",
    query: str = "",
    crs: str = "",
    **kwargs: Any,
) -> SpatVector:
    """
    Create a :class:`SpatVector`, like R ``terra::vect()``.

    * ``vect()`` — empty vector.
    * ``vect(str)`` — WKT literal or path to a vector file (via GDAL).
    * ``vect(SpatExtent)`` — rectangle as polygon (use **crs**).
    * ``list[str]`` — multiple WKT geometries.

    Extra GDAL arguments (``layer``, ``query``, …) match the C++ ``read`` call
    where applicable; see R ``terra::vect`` for full options (not all are wired yet).
    """
    del kwargs  # reserved for future parity

    if x is None:
        v = SpatVector()
        return messages(v, "vect")

    if isinstance(x, SpatExtent):
        crs_use = character_crs(crs, "vect") if crs else ""
        v = SpatVector(x, crs_use)
        return messages(v, "vect")

    if isinstance(x, str):
        s = x.strip()
        if _looks_like_wkt(s):
            v = SpatVector([s.replace("\n", "")])
            if crs:
                v.set_crs(character_crs(crs, "vect"))
            return messages(v, "vect")

        path = _normalize_path(s)
        v = SpatVector()
        ext: List[float] = []
        filt = SpatVector()
        opts: List[str] = []
        ok = v.read(path, layer, query, ext, filt, False, "", "", opts)
        if not ok:
            messages(v, "vect")
        if crs:
            v.set_crs(character_crs(crs, "vect"))
        return messages(v, "vect")

    if isinstance(x, (list, tuple)) and all(isinstance(i, str) for i in x):
        v = SpatVector([s.replace("\n", "") for s in x])
        if crs:
            v.set_crs(character_crs(crs, "vect"))
        return messages(v, "vect")

    raise TypeError(
        "vect: use None, str (path or WKT), list[str] (WKT), or SpatExtent"
    )
