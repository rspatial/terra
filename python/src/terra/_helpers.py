"""Shared helpers for the Python API (parallel to terra R helpers)."""

from __future__ import annotations

import warnings
from typing import Any, Union

__all__ = ["messages", "character_crs", "_getSpatDF", "_makeSpatDF"]


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


def _getSpatDF(sdf: Any) -> "Optional[pd.DataFrame]":
    """
    Convert a ``SpatDataFrame`` C++ object to a pandas DataFrame.

    Returns None if pandas is not available or the SpatDataFrame is empty.
    """
    try:
        import pandas as pd
    except ImportError:
        return None

    if sdf is None:
        return None

    # SpatDataFrame exposes .values() → Python dict {colname: list}
    if hasattr(sdf, "values"):
        d = sdf.values()
    elif isinstance(sdf, dict):
        d = sdf
    else:
        return None

    if not d:
        return pd.DataFrame()

    return pd.DataFrame(d)


def _makeSpatDF(df: "pd.DataFrame") -> Any:
    """
    Convert a pandas DataFrame to a ``SpatDataFrame`` C++ object.
    """
    import math
    from ._terra import SpatDataFrame

    sdf = SpatDataFrame()
    for col in df.columns:
        series = df[col]
        dtype = series.dtype
        kind = dtype.kind if hasattr(dtype, "kind") else "O"

        if kind == "f":                        # float64 / float32
            vals = []
            for v in series:
                try:
                    f = float(v)
                except (TypeError, ValueError):
                    f = float("nan")
                vals.append(f)
            sdf.add_column_double(vals, str(col))
        elif kind in ("i", "u"):               # int / uint
            # C++ long NA sentinel is the minimum long value
            _LONG_NA = -2147483648
            vals = []
            for v in series:
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    vals.append(_LONG_NA)
                else:
                    vals.append(int(v))
            sdf.add_column_long(vals, str(col))
        elif kind == "b":                      # bool
            vals = [int(bool(v)) if v is not None else 2 for v in series]
            sdf.add_column_bool(vals, str(col))
        else:                                  # string / object / other
            _STR_NA = "NA"
            vals = []
            for v in series:
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    vals.append(_STR_NA)
                else:
                    vals.append(str(v))
            sdf.add_column_string(vals, str(col))
    return sdf


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
