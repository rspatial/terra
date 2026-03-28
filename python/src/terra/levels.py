"""
levels.py — categorical raster levels (factors), color tables.
"""
from __future__ import annotations
from typing import Optional, Union, List, Dict, Any

from ._terra import SpatRaster, SpatOptions
from ._helpers import _getSpatDF, _makeSpatDF, messages
from .names import set_names_inplace


def _opt() -> SpatOptions:
    return SpatOptions()


def _layer_idx(x: SpatRaster, layer: Union[int, str]) -> int:
    """0-based layer index (matches C++ ``source`` layer indexing)."""
    if isinstance(layer, str):
        return list(x.names).index(layer)
    return int(layer)


# ---------------------------------------------------------------------------
# is_factor / as_factor
# ---------------------------------------------------------------------------

def is_factor(x: SpatRaster) -> List[bool]:
    """Return per-layer flags indicating whether the layer is categorical."""
    return list(x.hasCategories())


def as_factor(x: SpatRaster) -> SpatRaster:
    """Convert *x* to a categorical raster (in-place on a copy)."""
    x = x.deepcopy() if hasattr(x, 'deepcopy') else x
    x.pntr = x.makeCategorical(-1, _opt())
    return messages(x, "as_factor")


# ---------------------------------------------------------------------------
# levels / set levels
# ---------------------------------------------------------------------------

def levels(x: SpatRaster, layer: Optional[Union[int, str]] = None) -> List[Any]:
    """
    Return the categorical lookup table(s) for *x*.

    Parameters
    ----------
    x : SpatRaster
    layer : int or str, optional
        If given, return only the table for that layer.

    Returns
    -------
    list of pandas.DataFrames (or None for non-factor layers).
    """
    cats = x.getCategories()
    result = []
    for cat in cats:
        df = _getSpatDF(cat.df)
        if df is None or len(df.columns) == 0:
            result.append(None)
        else:
            idx = cat.index + 1
            cols = [df.columns[0], df.columns[min(idx, len(df.columns) - 1)]]
            result.append(df[cols])
    if layer is not None:
        if isinstance(layer, str):
            nms = list(x.names)
            layer = nms.index(layer)
        return result[layer]
    return result


def set_levels(x: SpatRaster, value: Union[List, Any], active: int = 1) -> SpatRaster:
    """
    Set the categorical lookup table for *x*.

    Parameters
    ----------
    x : SpatRaster
    value : DataFrame or list of DataFrames
        Lookup table(s) with columns [value_id, label, ...].
    active : int
        Which column (1-based) is the active label column.

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    if not isinstance(value, list):
        value = [value]
    for i, v in enumerate(value):
        if v is None:
            xc.removeCategories(i)
        else:
            set_cats(xc, i, v, active)
    return xc


# ---------------------------------------------------------------------------
# cats / set_cats / categories
# ---------------------------------------------------------------------------

def cats(x: SpatRaster, layer: Optional[Union[int, str]] = None) -> List[Any]:
    """
    Return the full category data frame(s) for *x*.

    Parameters
    ----------
    x : SpatRaster
    layer : int or str, optional
        If given, return only for that layer.

    Returns
    -------
    list of pandas.DataFrames (or None).
    """
    raw = x.getCategories()
    result = []
    for cat in raw:
        if cat.df.nrow == 0:
            result.append(None)
        else:
            result.append(_getSpatDF(cat.df))
    if layer is not None:
        if isinstance(layer, str):
            nms = list(x.names)
            layer = nms.index(layer)
        return result[layer]
    return result


def set_cats(
    x: SpatRaster,
    layer: Union[int, str] = 0,
    value: Any = None,
    active: int = 1,
) -> bool:
    """
    Set the category data frame for a single layer of *x* (in-place).

    Parameters
    ----------
    x : SpatRaster
    layer : int or str
        Layer index (0-based) or name.
    value : pandas.DataFrame
        Lookup table with at least two columns: integer IDs and labels.
    active : int
        Which label column (1-based offset past the ID column) is active.

    Returns
    -------
    bool
    """
    layer_idx = _layer_idx(x, layer)

    if value is None:
        x.removeCategories(layer_idx)
        return True

    import pandas as pd
    if not isinstance(value, pd.DataFrame):
        raise TypeError("value must be a pandas DataFrame")
    if len(value.columns) < 2:
        raise ValueError("value must have at least two columns")
    # Ensure first col is int
    value = value.copy()
    value.iloc[:, 0] = value.iloc[:, 0].astype(int)
    sdf = _makeSpatDF(value)
    idx = max(1, min(len(value.columns) - 1, active))
    ok = x.setCategories(layer_idx, sdf, idx)
    return bool(ok)


def categories(
    x: SpatRaster,
    layer: Union[int, str] = 0,
    value: Any = None,
    active: int = 1,
) -> SpatRaster:
    """
    Return a copy of *x* with updated category table for *layer*.

    Parameters
    ----------
    x : SpatRaster
    layer : int or str
    value : pandas.DataFrame
    active : int

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    set_cats(xc, layer, value, active)
    return xc


# ---------------------------------------------------------------------------
# active_cat
# ---------------------------------------------------------------------------

def active_cat(x: SpatRaster, layer: Union[int, str] = 0) -> int:
    """
    Return the index of the active label column for a categorical layer.

    Parameters
    ----------
    x : SpatRaster
    layer : int or str
        Layer (0-based index or name).

    Returns
    -------
    int
        Column index within the category table (C++ convention).
    """
    return int(x.getCatIndex(_layer_idx(x, layer)))


def set_active_cat(
    x: SpatRaster,
    value: Union[int, str],
    layer: Union[int, str] = 0,
) -> SpatRaster:
    """
    Return a copy of *x* with a different active category column.

    Parameters
    ----------
    x : SpatRaster
    value : int or str
        Column index or column name to activate.
    layer : int or str
        Layer (0-based index or name).

    Returns
    -------
    SpatRaster
    """
    layer_idx = _layer_idx(x, layer)

    if isinstance(value, str):
        ct = cats(x, layer_idx)
        if ct is None:
            raise ValueError("layer is not categorical")
        col_idx = list(ct.columns).index(value)
    else:
        col_idx = int(value)

    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    if not xc.setCatIndex(layer_idx, col_idx):
        raise RuntimeError("invalid category index")
    return xc


# ---------------------------------------------------------------------------
# add_cats / concats / catalyze / drop_levels
# ---------------------------------------------------------------------------

def add_cats(
    x: SpatRaster,
    value: Any,
    merge: bool = False,
    layer: Union[int, str] = 0,
) -> SpatRaster:
    """
    Add columns to the category table of a layer.

    Parameters
    ----------
    x : SpatRaster
    value : pandas.DataFrame
        New column(s) to add; first column must be the ID.
    merge : bool
        If True, merge by ID rather than binding columns.
    layer : int or str

    Returns
    -------
    SpatRaster
    """
    import pandas as pd

    li = _layer_idx(x, layer)
    cts = cats(x, li)
    if cts is None:
        raise ValueError("layer has no categories to add to")
    nact = len(cts.columns)
    if merge:
        if len(value.columns) < 2:
            raise ValueError("value must have at least two columns when merge=True")
        cts = pd.merge(cts, value, on=cts.columns[0], how="left")
        cts = cts.sort_values(cts.columns[0]).reset_index(drop=True)
    else:
        if len(cts) != len(value):
            raise ValueError("number of categories does not match")
        cts = pd.concat([cts, value.reset_index(drop=True)], axis=1)
    return categories(x, layer=li, value=cts, active=nact)


def drop_levels(x: SpatRaster) -> SpatRaster:
    """Drop unused levels from all categorical layers."""
    xc = x.droplevels()
    return messages(xc, "drop_levels")


def concats(x: SpatRaster, y: SpatRaster) -> SpatRaster:
    """Combine the category tables of *x* and *y*."""
    opt = _opt()
    xc = x.combineCats(y, opt)
    return messages(xc, "concats")


def catalyze(x: SpatRaster) -> SpatRaster:
    """Convert categorical column values to separate numeric layers."""
    ct = cats(x)
    if not any(c is not None for c in ct):
        return x
    from_vals = []
    to_vals_list = []
    nms = None
    for i, c in enumerate(ct):
        if c is not None and len(c) > 0:
            nms = list(c.columns[1:])
            from_vals = c.iloc[:, 0].tolist()
            for j in range(1, len(c.columns)):
                col = c.iloc[:, j]
                try:
                    to_vals_list.append([float(v) for v in col])
                except (ValueError, TypeError):
                    import pandas as pd
                    to_vals_list.append([float(i) for i in pd.Categorical(col).codes + 1])
            break
    if not from_vals:
        return x
    opt = _opt()
    xc = x.lookup_catalyze(from_vals, to_vals_list, opt)
    result = messages(xc, "catalyze")
    if nms:
        set_names_inplace(result, [str(n) for n in nms])
    return result


# ---------------------------------------------------------------------------
# color tables
# ---------------------------------------------------------------------------

def has_colors(x: SpatRaster) -> List[bool]:
    """Return per-layer flags indicating whether a color table exists."""
    return list(x.hasColors())


def coltab(x: SpatRaster) -> List[Any]:
    """
    Return the color table(s) for *x*.

    Returns
    -------
    list of pandas.DataFrames (or None for layers without color tables).
    """
    has = list(x.hasColors())
    if not any(has):
        return [None] * len(has)
    raw = x.getColors()
    result = []
    for i, df_raw in enumerate(raw):
        if has[i]:
            result.append(_getSpatDF(df_raw))
        else:
            result.append(None)
    return result


def set_coltab(
    x: SpatRaster,
    value: Any,
    layer: Union[int, str] = 0,
) -> SpatRaster:
    """
    Set a color table for a layer of *x*.

    Parameters
    ----------
    x : SpatRaster
    value : list of str (hex colors), array, or DataFrame with RGBA columns.
    layer : int or str

    Returns
    -------
    SpatRaster
    """
    import pandas as pd
    import numpy as np

    layer_idx = _layer_idx(x, layer)

    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x

    if value is None:
        xc.removeColors(layer_idx)
        return xc

    if isinstance(value, (list, np.ndarray)) and not isinstance(value, pd.DataFrame):
        value_list = list(value)
        if len(value_list) > 0 and isinstance(value_list[0], str):
            # Parse hex strings
            rows = []
            for idx, hex_color in enumerate(value_list):
                hex_color = hex_color.lstrip('#')
                if len(hex_color) == 6:
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    a = 255
                elif len(hex_color) == 8:
                    r, g, b, a = (int(hex_color[i:i+2], 16) for i in range(0, 8, 2))
                else:
                    r, g, b, a = 128, 128, 128, 255
                rows.append([idx, r, g, b, a])
            value = pd.DataFrame(rows, columns=["value", "red", "green", "blue", "alpha"])

    if isinstance(value, pd.DataFrame):
        value = value.copy()
        value.iloc[:, 0] = value.iloc[:, 0].astype(int)
        for col in value.columns[1:]:
            value[col] = value[col].clip(0, 255).astype(int)
        value = value.fillna(255)
        sdf = _makeSpatDF(value)
        if not xc.setColors(layer_idx, sdf):
            raise RuntimeError("could not set color table")
    return xc
