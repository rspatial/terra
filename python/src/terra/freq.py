"""
freq — frequency table of raster values (R ``terra::freq`` for ``SpatRaster``).
"""
from __future__ import annotations

import warnings
from typing import Any, List, Optional, Union

import numpy as np

from ._terra import SpatOptions, SpatRaster
from ._helpers import messages
from .levels import active_cat, cats, is_factor


def _freq_apply_factor_labels(x: SpatRaster, df: "pd.DataFrame") -> "pd.DataFrame":
    """Replace coded ``value`` with category labels for factor layers (R ``cats`` loop)."""
    ff = is_factor(x)
    if not any(ff):
        return df
    out = df.copy()
    cgs = cats(x)
    for f in range(1, x.nlyr() + 1):
        if not ff[f - 1]:
            continue
        cg = cgs[f - 1]
        if cg is None or len(cg.columns) < 2:
            continue
        mask = out["layer"] == f
        if not mask.any():
            continue
        ac = active_cat(x, f)
        code_col = cg.columns[0]
        lab_col = cg.columns[min(ac, len(cg.columns) - 1)]
        for idx in out.loc[mask].index:
            v = out.at[idx, "value"]
            row = cg[cg[code_col] == v]
            if len(row) == 0 and isinstance(v, (int, float, np.floating)):
                fv = float(v)
                if np.isnan(fv):
                    row = cg[np.isnan(cg[code_col].astype(float))]
                else:
                    row = cg[np.isclose(cg[code_col].astype(float), fv)]
            if len(row):
                out.at[idx, "value"] = row[lab_col].iloc[0]
    return out


def _freq_cpp_to_dataframe(
    x: SpatRaster,
    raw: List[List[float]],
    *,
    bylayer: bool,
) -> "pd.DataFrame":
    """Turn C++ ``freq()`` output into a long data frame."""
    import pandas as pd

    rows: List[tuple] = []

    if bylayer:
        for lyr in range(len(raw)):
            vec = raw[lyr]
            if not vec:
                continue
            arr = np.asarray(vec, dtype=float)
            if arr.size == 0 or np.all(np.isnan(arr)):
                continue
            n = len(arr) // 2
            vals = arr[:n]
            cnts = arr[n:]
            for v, c in zip(vals, cnts):
                rows.append((lyr + 1, float(v), float(c)))
    else:
        vec = raw[0] if raw else []
        arr = np.asarray(vec, dtype=float) if len(vec) else np.array([])
        if arr.size:
            n = len(arr) // 2
            vals = arr[:n]
            cnts = arr[n:]
            for v, c in zip(vals, cnts):
                rows.append((float(v), float(c)))

    if not rows:
        if bylayer:
            return pd.DataFrame(columns=["layer", "value", "count"])
        return pd.DataFrame(columns=["value", "count"])

    if bylayer:
        return pd.DataFrame(rows, columns=["layer", "value", "count"])
    return pd.DataFrame(rows, columns=["value", "count"])


def freq(
    x: SpatRaster,
    digits: Optional[int] = 0,
    value: Any = None,
    bylayer: bool = True,
    usenames: bool = False,
    zones: Any = None,
    wide: bool = False,
    touches: bool = False,
) -> Any:
    """
    Frequency table of cell values in a :class:`SpatRaster`.

    Mirrors R ``terra::freq`` for ``SpatRaster`` (``y`` missing).

    Parameters
    ----------
    x : SpatRaster
    digits : int or None
        ``None`` = no rounding (R ``digits=NA``).  An ``int`` (default ``0``)
        turns on rounding with that many decimal places.
    value : scalar, optional
        If set, return only the total count of this value across cells
        (per layer if ``bylayer``).  Use ``float('nan')`` to count NAs.
    bylayer : bool
    usenames : bool
        Replace layer indices with layer names.
    zones : optional
        Not implemented in Python (raises ``NotImplementedError``).
    wide : bool
        Reshape counts to wide format (one column per distinct value).
    touches : bool
        Reserved for future ``zones`` support (vector).

    Returns
    -------
    pandas.DataFrame
    """
    import pandas as pd

    if zones is not None:
        raise NotImplementedError(
            "freq(..., zones=...) is not implemented in Python yet."
        )

    if not x.hasValues:
        warnings.warn("SpatRaster x has no cell values", UserWarning, stacklevel=2)
        return pd.DataFrame(columns=["layer", "value", "count"])

    opt = SpatOptions()
    nms = list(x.names)

    # --- count a single value -------------------------------------------------
    if value is not None:
        if isinstance(value, str):
            ff = is_factor(x)
            if not any(ff):
                raise ValueError(
                    "a character value is only meaningful for categorical rasters"
                )
            fidx = next(i for i, b in enumerate(ff) if b)
            sub = messages(x.subset([fidx], opt), "subset")
            df0 = freq(sub, digits=digits, value=None, bylayer=True, usenames=False)
            df0 = df0[df0["layer"] == 1].copy()
            out = df0[df0["value"].astype(str) == value]
            if usenames:
                out = out.copy()
                out["layer"] = nms[fidx]
            return out.reset_index(drop=True)

        use_round = digits is not None
        d = 0 if digits is None else int(digits)
        vfloat = float(value)
        cnts = x.count(vfloat, bylayer, use_round, d, opt)
        cnts = [int(c) for c in cnts]
        if np.isnan(vfloat):
            v_out = float("nan")
        elif use_round:
            v_out = float(round(vfloat, d)) if d > 0 else float(round(vfloat))
        else:
            v_out = vfloat

        if bylayer:
            df = pd.DataFrame(
                {
                    "layer": np.arange(1, x.nlyr() + 1),
                    "value": v_out,
                    "count": cnts,
                }
            )
        else:
            df = pd.DataFrame({"value": [v_out], "count": [int(sum(cnts))]})
        if usenames and bylayer:
            df = df.copy()
            df["layer"] = [nms[i - 1] for i in df["layer"]]
        return df

    # --- full table -----------------------------------------------------------
    use_round = digits is not None
    d = 0 if digits is None else int(digits)
    raw = x.freq(bylayer, use_round, d, opt)

    df = _freq_cpp_to_dataframe(x, raw, bylayer=bylayer)

    if bylayer and len(df) > 0:
        df = _freq_apply_factor_labels(x, df)
    elif (not bylayer) and len(df) > 0 and any(is_factor(x)):
        # Pooled counts: cannot map per-layer category codes; leave numeric.
        pass

    if usenames and bylayer and len(df) > 0:
        df = df.copy()
        df["layer"] = df["layer"].astype(int).map(lambda i: nms[i - 1])

    if wide and len(df) > 0:
        df = df.copy()
        df["count"] = df["count"].fillna(0)
        if bylayer:
            pv = df.pivot(index="layer", columns="value", values="count").fillna(0)
            pv.columns.name = None
            return pv.reset_index()
        pv = df.pivot(columns="value", values="count").fillna(0)
        return pv

    df = df[df["count"].notna()]
    return df.reset_index(drop=True)
