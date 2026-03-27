"""
extract.py — extract raster values at points, lines, or polygons.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union

from ._terra import SpatRaster, SpatVector
from ._helpers import messages, spatoptions


def _flat_extract_to_dataframe(
    flat: List[float],
    *,
    geo: str,
    nl: int,
    nc: int,
    cn: List[str],
    cells: bool,
) -> "pd.DataFrame":
    """Reshape ``extractVectorFlat`` output like R ``extract.R`` (matrix + colnames)."""
    import numpy as np
    import pandas as pd

    e = np.asarray(flat, dtype=float)
    if e.size == 0:
        return pd.DataFrame(columns=cn)

    if geo == "points":
        if nc == nl:
            nrow = e.size // nc
            mat = e.reshape((nc, nrow), order="F").T
        else:
            nrow = e.size // nc
            mat = e.reshape((nrow, nc), order="C")
        ids = np.arange(1, nrow + 1, dtype=float)
        arr = np.column_stack([ids, mat])
    else:
        ncol = nc + 1
        nrow = e.size // ncol
        arr = e.reshape((nrow, ncol), order="C")

    df = pd.DataFrame(arr, columns=list(cn))
    if cells and "cell" in df.columns:
        df["cell"] = df["cell"] + 1.0
    return df


def extract(
    x: SpatRaster,
    y: Union[SpatVector, "np.ndarray", "pd.DataFrame", List],
    *,
    method: str = "simple",
    layer: Optional[Union[int, str, "SpatRaster"]] = None,
    fun: Optional[Union[str, Callable]] = None,
    na_rm: bool = True,
    bind: bool = False,
    raw: bool = False,
    cells: bool = False,
    xy: bool = False,
    ID: bool = True,
    weights: bool = False,
    exact: bool = False,
    touches: Optional[bool] = None,
    small: bool = True,
    factors: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> "pd.DataFrame":
    """
    Extract raster values at spatial locations.

    Parameters
    ----------
    x : SpatRaster
        Raster with values to extract.
    y : SpatVector, array, or list
        Locations.  Points, lines, or polygons.  Can also be a matrix or
        list of ``[x_coord, y_coord]`` pairs.
    method : str
        ``"simple"`` (cell-centre, default) or ``"bilinear"``
        (4-cell bilinear interpolation; points only).
    layer : int, str, or SpatRaster, optional
        Variable-layer extraction: start layer index per feature (1-based),
        or a SpatRaster whose values give per-cell start indices.
    fun : str or callable, optional
        Aggregation function for polygons/lines.  If None, all intersecting
        values are returned.
    na_rm : bool
        Ignore NA values when aggregating.
    bind : bool
        Append extracted values to the attribute table of *y*.
    raw : bool
        Return numeric matrix without conversion of factor levels.
    cells : bool
        Include cell numbers in the output.
    xy : bool
        Include x/y coordinates in the output.
    ID : bool
        Include feature IDs in the output.
    weights : bool
        Include fractional cell coverage weights (polygons only).
    exact : bool
        Compute exact fractional weights (slower).
    touches : bool, optional
        Include cells that touch the polygon boundary.
    small : bool
        Include cells where the polygon is smaller than a cell.
    factors : bool
        Return factor levels as strings rather than integers.
    filename : str
    overwrite : bool

    Returns
    -------
    pandas.DataFrame
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for extract()")
    from ._helpers import _getSpatDF

    # Coerce y to SpatVector if needed
    if not isinstance(y, SpatVector):
        from .vect import vect
        y = vect(y)

    opt = spatoptions(filename, overwrite)

    if layer is not None:
        raise NotImplementedError("extract(..., layer=...) is not implemented in the Python API yet")

    if fun is not None:
        fun_str = fun if isinstance(fun, str) else getattr(fun, "__name__", "mean")
        res_sv = x.extractByValues(y, fun_str, na_rm, opt)
        messages(x, "extract")
        df = _getSpatDF(res_sv)
        if df is None:
            return pd.DataFrame()
        if ID:
            df.insert(0, "ID", range(1, len(df) + 1))
        return df

    if touches is None:
        touches = True

    geo = y.type()
    flat = x.extractVectorFlat(y, [""], False, touches, small, method, cells, xy, weights, exact, opt)
    messages(x, "extract")

    nl = x.nlyr()
    nc = nl
    names = list(x.names)
    cn: List[str] = ["ID"] + names
    if cells:
        cn.append("cell")
        nc += 1
    if xy:
        cn.extend(["x", "y"])
        nc += 2
    if weights:
        cn.append("weight")
        nc += 1
    elif exact:
        cn.append("fraction")
        nc += 1

    df = _flat_extract_to_dataframe(flat, geo=geo, nl=nl, nc=nc, cn=cn, cells=cells)

    if not ID and "ID" in df.columns:
        df = df.drop(columns=["ID"])

    if factors:
        cats = x.getCategories()
        for i, cat in enumerate(cats):
            lnm = list(x.names)[i]
            ct_df = _getSpatDF(cat.df) if cat.df.nrow > 0 else None
            if ct_df is not None and lnm in df.columns:
                idx_col = ct_df.columns[0]
                lbl_col = ct_df.columns[cat.index + 1] if cat.index + 1 < len(ct_df.columns) else ct_df.columns[-1]
                mapping = dict(zip(ct_df[idx_col], ct_df[lbl_col]))
                df[lnm] = df[lnm].map(mapping)

    if bind:
        y_df = _getSpatDF(y.df)
        if y_df is not None and len(y_df) == len(df):
            df = pd.concat([y_df.reset_index(drop=True), df.reset_index(drop=True)], axis=1)

    return df


def extract_xy(
    x: SpatRaster,
    xy: Union["np.ndarray", List],
    method: str = "simple",
) -> "pd.DataFrame":
    """
    Extract raster values at x/y coordinate pairs.

    Parameters
    ----------
    x : SpatRaster
    xy : array-like, shape (n, 2)
        Coordinate pairs (x_col, y_col).
    method : str
        ``"simple"`` or ``"bilinear"``.

    Returns
    -------
    pandas.DataFrame with columns matching layer names.
    """
    import numpy as np
    import pandas as pd
    from .vect import vect

    xy_arr = np.asarray(xy, dtype=float)
    if xy_arr.ndim == 1:
        xy_arr = xy_arr.reshape(1, -1)
    if xy_arr.shape[1] < 2:
        raise ValueError("xy must have at least 2 columns")

    pts = vect(xy_arr[:, :2], type="points", crs=x.crs)
    return extract(x, pts, method=method)
