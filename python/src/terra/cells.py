"""
cells.py — cell index operations for SpatRaster.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatExtent, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# cells
# ---------------------------------------------------------------------------

def cells(
    x: SpatRaster,
    y: Optional[Union[float, List[float], SpatVector, SpatExtent]] = None,
    *,
    method: str = "simple",
    weights: bool = False,
    exact: bool = False,
    touches: Optional[bool] = None,
    small: bool = True,
    pairs: bool = False,
) -> Union[np.ndarray, List[np.ndarray]]:
    """
    Return cell numbers for non-NA cells or cells matching given criteria.

    Parameters
    ----------
    x : SpatRaster
    y : float, list, SpatVector, or SpatExtent, optional
        - If absent: return non-NA cell indices (1-based).
        - If a scalar or list: find cells whose values match *y*.
        - If a SpatVector: find cells overlapping *y*.
        - If a SpatExtent: find cells within *y*.
    method : str
        ``"simple"`` or ``"bilinear"`` (for SpatVector point lookups).
    weights : bool
        Include fractional coverage weights (SpatVector polygons only).
    exact : bool
        Compute exact fractional weights.
    touches : bool, optional
        Include cells touching the polygon boundary.
    small : bool
        Include small polygons.
    pairs : bool
        If *y* is numeric, return ``(cell, value)`` pairs.

    Returns
    -------
    numpy.ndarray (1-D cell indices, 1-based), or
    list of arrays (one per layer when *y* is numeric), or
    numpy.ndarray with columns [ID, cell, …] when *y* is a SpatVector.
    """
    opt = _opt()

    if y is None:
        raw = x.cells_notna_novalues(opt)
        return np.array(raw, dtype=int) + 1

    if isinstance(y, SpatExtent):
        raw = x.extCells(y.pntr) if hasattr(y, 'pntr') else x.extCells(y)
        return np.array(raw, dtype=int) + 1

    if isinstance(y, SpatVector):
        if touches is None:
            from .vect import is_lines
            touches = is_lines(y)
        raw = x.vectCells(y, touches, small, method, weights, exact, opt)
        if y.geomtype() == "points":
            if method == "bilinear":
                m = np.array(raw, dtype=float).reshape(y.nrow(), -1)
                m[:, 0:4] = m[:, 0:4] + 1
                ids = np.arange(1, y.nrow() + 1).reshape(-1, 1)
                m = np.hstack([ids, m])
                return m
            else:
                m = np.array(raw, dtype=float).reshape(y.nrow(), -1)
                m[:, 0] = m[:, 0] + 1
                ids = np.arange(1, y.nrow() + 1).reshape(-1, 1)
                m = np.hstack([ids, m])
                return m
        else:
            ncols = 3 if (weights or exact) else 2
            m = np.array(raw, dtype=float).reshape(-1, ncols)
            m[:, :2] = m[:, :2] + 1
            return m

    # Numeric value matching
    if isinstance(y, (int, float)):
        y = [float(y)]
    else:
        y = [float(v) for v in y]
    raw_list = x.is_in_cells(y, pairs, opt)
    result = messages(x, "cells")
    if pairs:
        out = []
        for arr in raw_list:
            m = np.array(arr, dtype=float).reshape(-1, 2)
            m[:, 0] = m[:, 0] + 1
            out.append(m)
        return out
    else:
        return [np.array(arr, dtype=int) + 1 for arr in raw_list]


# ---------------------------------------------------------------------------
# Row / col / cell conversions
# ---------------------------------------------------------------------------

def row_from_y(x: SpatRaster, y: Union[float, List[float]]) -> np.ndarray:
    """Return row indices (1-based) for y-coordinates."""
    if isinstance(y, (int, float)):
        y = [float(y)]
    return np.array(x.rowFromY(y), dtype=int) + 1


def col_from_x(x: SpatRaster, xcoord: Union[float, List[float]]) -> np.ndarray:
    """Return column indices (1-based) for x-coordinates."""
    if isinstance(xcoord, (int, float)):
        xcoord = [float(xcoord)]
    return np.array(x.colFromX(xcoord), dtype=int) + 1


def cell_from_xy(
    x: SpatRaster,
    xy: Union[List, np.ndarray],
) -> np.ndarray:
    """
    Return cell numbers for x/y coordinate pairs.

    Parameters
    ----------
    x : SpatRaster
    xy : array-like, shape (n, 2)
        Columns: [x_coord, y_coord].

    Returns
    -------
    numpy.ndarray, dtype float64
        1-based cell numbers.  Invalid coordinates (outside extent or NaN
        inputs) are ``nan``, matching R's ``NA`` for cell index.
    """
    xy = np.asarray(xy, dtype=float)
    if xy.ndim == 1:
        xy = xy.reshape(1, -1)
    raw = x.cellFromXY(xy[:, 0].tolist(), xy[:, 1].tolist(), float("nan"))
    return np.array(raw, dtype=float) + 1.0


def cell_from_row_col(
    x: SpatRaster,
    row: Union[int, List[int]],
    col: Union[int, List[int]],
) -> np.ndarray:
    """Return cell numbers for row/column index pairs (1-based)."""
    if isinstance(row, int):
        row = [row]
    if isinstance(col, int):
        col = [col]
    row0 = [r - 1 for r in row]
    col0 = [c - 1 for c in col]
    raw = x.cellFromRowCol(row0, col0)
    return np.array(raw, dtype=int) + 1


def xy_from_cell(x: SpatRaster, cell: Union[int, List[int]]) -> np.ndarray:
    """
    Return x/y coordinates for cell numbers.

    Parameters
    ----------
    x : SpatRaster
    cell : int or list of int (1-based)

    Returns
    -------
    numpy.ndarray, shape (n, 2), columns [x, y].
    """
    if isinstance(cell, int):
        cell = [cell]
    cell0 = [c - 1 for c in cell]
    coords = x.xyFromCell(cell0)
    return np.array(coords, dtype=float).reshape(-1, 2)


def row_col_from_cell(
    x: SpatRaster,
    cell: Union[int, List[int]],
) -> np.ndarray:
    """
    Return row and column indices for cell numbers.

    Parameters
    ----------
    x : SpatRaster
    cell : int or list (1-based)

    Returns
    -------
    numpy.ndarray, shape (n, 2), columns [row, col] (1-based).
    """
    if isinstance(cell, int):
        cell = [cell]
    cell0 = [c - 1 for c in cell]
    rc = x.rowColFromCell(cell0)
    arr = np.array(rc, dtype=int).reshape(-1, 2) + 1
    return arr
