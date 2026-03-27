"""
crosstab — cross-tabulate raster layers (R ``terra::crosstab``).
"""
from __future__ import annotations

import re
from typing import Any, List

import numpy as np

from ._terra import SpatRaster, SpatOptions
from .levels import is_factor, levels as rast_levels


def _make_names_unique(names: List[str]) -> List[str]:
    """Approximate R ``make.names(..., unique=TRUE)`` for column / index labels."""
    out: List[str] = []
    seen: dict[str, bool] = {}
    for raw in names:
        s = re.sub(r"[^A-Za-z0-9_.]", ".", str(raw)) or "X"
        if s[0].isdigit():
            s = "X" + s
        base = s
        cand = base
        i = 0
        while cand in seen:
            i += 1
            cand = f"{base}.{i}"
        seen[cand] = True
        out.append(cand)
    return out


def crosstab(
    x: SpatRaster,
    digits: int = 0,
    long: bool = False,
    useNA: bool = False,
) -> Any:
    """
    Cross-tabulate the layers of a :class:`SpatRaster` (contingency table).

    Parameters
    ----------
    x : SpatRaster
        Must have at least two layers.
    digits : int
        Decimal places for rounding values before counting (``digits < 0`` skips
        rounding in the C++ implementation).
    long : bool
        If True, return a long data frame with counts in column ``"n"``.
        If False, return a wide table (two layers → 2-D DataFrame; more layers →
        Series with :class:`pandas.MultiIndex`).
    useNA : bool
        If False (default), cells where any layer is NA are dropped (passed to C++
        as ``narm=True``). If True, NA is treated as a category.

    Returns
    -------
    pandas.DataFrame or pandas.Series
        Wide format is a contingency table; long format has one row per unique
        combination with positive count.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError("pandas is required for crosstab()") from e

    nl = x.nlyr()
    if nl < 2:
        raise ValueError("crosstab needs at least 2 layers")

    nms = list(x.names)
    opt = SpatOptions()
    res = x.crosstab(int(digits), not useNA, opt)

    if len(res) == 0:
        return pd.DataFrame(columns=nms + ["Freq"])

    arr = np.asarray(res, dtype=float).reshape(-1, nl + 1)
    df = pd.DataFrame(arr, columns=nms + ["Freq"])

    # Replace factor codes with labels (R is.factor / levels)
    ff = is_factor(x)
    if any(ff):
        lvl_list = rast_levels(x)
        for i, is_f in enumerate(ff):
            if not is_f:
                continue
            v = lvl_list[i]
            if v is None or len(v.columns) < 2:
                continue
            code_col, lab_col = v.columns[0], v.columns[1]
            tab = v.set_index(code_col)[lab_col]
            df[nms[i]] = df[nms[i]].map(tab)

    if long:
        df = df.loc[df["Freq"] > 0].copy()
        df = df.rename(columns={"Freq": "n"})
        return df.reset_index(drop=True)

    safe = _make_names_unique(nms)
    df.columns = safe + ["Freq"]

    if nl == 2:
        out = df.pivot(index=safe[0], columns=safe[1], values="Freq")
        if not useNA:
            out = out.fillna(0)
        out = out.sort_index().sort_index(axis=1)
        return out

    # Multi-way table: one row per combination, MultiIndex over layer columns
    ser = df.set_index(safe)["Freq"]
    ser = ser.sort_index()
    return ser
