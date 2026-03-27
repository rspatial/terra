"""CRS helpers — parallel to R ``terra::crs()``."""

from __future__ import annotations

from typing import Any, Optional

from ._helpers import character_crs, messages

__all__ = ["crs"]


def crs(x: Any, value: Optional[str] = None, *, proj4: bool = False) -> Any:
    """
    Get or set the coordinate reference system, like R ``terra::crs()``.

    * ``crs(x)`` — return CRS string (WKT by default; use ``proj4=True`` for PROJ.4).
    * ``crs(x, value)`` — set CRS from a string or another object; returns *x*.
    """
    if value is not None:
        s = character_crs(value, "crs")
        x.set_crs(s)
        return messages(x, "crs")

    kind = "proj4" if proj4 else "wkt"
    return x.get_crs(kind)
