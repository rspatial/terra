"""Shared helpers for the Python API (parallel to terra R helpers)."""

from __future__ import annotations

import warnings
from typing import Any, Union

__all__ = ["messages", "character_crs"]


def messages(obj: Any, caller: str = "") -> Any:
    """
    Surface C++ warnings and errors stored on *obj*, like R's ``messages()``.

    Issues Python :class:`warnings.warn` for warnings and raises
    :exc:`RuntimeError` on error. Returns *obj* unchanged if there is no error.
    """
    prefix = f"[{caller}] " if caller else ""

    if hasattr(obj, "has_warning") and obj.has_warning():
        wmsg = " ".join(obj.getWarnings())
        warnings.warn(f"{prefix}{wmsg}", stacklevel=2)

    if hasattr(obj, "has_error") and obj.has_error():
        raise RuntimeError(f"{prefix}{obj.getError()}")

    return obj


def character_crs(x: Union[str, Any], _caller: str = "") -> str:
    """
    Normalise a CRS string similarly to R's ``character_crs`` (subset).
    """
    if x is None:
        return ""
    if isinstance(x, str):
        s = x
    else:
        # Another object with crs — use get_crs if present
        if hasattr(x, "get_crs"):
            return character_crs(x.get_crs("wkt"), _caller)
        s = str(x)

    if not s or s.lower() == "nan":
        return ""
    t = s.strip().lower()
    if t == "local":
        return (
            'LOCAL_CS["Cartesian (Meter)", LOCAL_DATUM["Local Datum",0], '
            'UNIT["Meter",1.0], AXIS["X",EAST], AXIS["Y",NORTH]]'
        )
    if t == "lonlat":
        return "+proj=longlat"
    return s
