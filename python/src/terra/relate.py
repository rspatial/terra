"""
relate.py — spatial relationships between SpatVector geometries.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatExtent, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


_DE9IM_RELATIONS = {
    "intersects", "touches", "crosses", "overlaps",
    "within", "contains", "covers", "coveredBy",
    "disjoint", "equals",
}


def _to_vect(x) -> SpatVector:
    """Convert SpatRaster or SpatExtent to a SpatVector of polygons."""
    from .vect import vect as make_vect
    if isinstance(x, SpatVector):
        return x
    if isinstance(x, SpatRaster):
        return x.as_polygons(True, False, True, False, _opt())
    if isinstance(x, SpatExtent):
        return x.as_polygons()
    raise TypeError(f"Cannot convert {type(x)} to SpatVector")


def is_related(
    x: Union[SpatVector, SpatRaster, SpatExtent],
    y: Union[SpatVector, SpatRaster, SpatExtent],
    relation: str,
) -> np.ndarray:
    """
    Test a binary spatial relation between *x* and *y*.

    Parameters
    ----------
    x, y : SpatVector, SpatRaster, or SpatExtent
    relation : str
        One of: ``"intersects"``, ``"touches"``, ``"crosses"``,
        ``"overlaps"``, ``"within"``, ``"contains"``, ``"covers"``,
        ``"coveredBy"``, ``"disjoint"``, ``"equals"``.

    Returns
    -------
    numpy.ndarray of bool, one value per feature in *x*.
    """
    xv = _to_vect(x)
    yv = _to_vect(y)
    out = xv.is_related(yv, relation)
    messages(xv, "is_related")
    return np.array(out, dtype=bool)


def relate(
    x: SpatVector,
    y: SpatVector,
    relation: str,
    *,
    pairs: bool = False,
    na_rm: bool = True,
) -> Union[np.ndarray, "pd.DataFrame"]:
    """
    Compute a spatial relation matrix between *x* and *y*.

    Parameters
    ----------
    x : SpatVector
    y : SpatVector
    relation : str
        DE-9IM relation name or pattern.
    pairs : bool
        If True, return a 2-column array of feature index pairs (1-based)
        instead of a boolean matrix.
    na_rm : bool
        Remove NA results when *pairs=True*.

    Returns
    -------
    numpy.ndarray, shape (nrow(x), nrow(y)) of bool, or
    numpy.ndarray, shape (n, 2) of [id.x, id.y] (1-based) if pairs=True.
    """
    out = x.related_between(y, relation, na_rm)
    messages(x, "relate")
    if pairs:
        if len(out[0]) == 0:
            return np.empty((0, 2), dtype=int)
        m = np.column_stack([np.array(out[0], dtype=int) + 1,
                             np.array(out[1], dtype=int) + 1])
        return m
    else:
        m = np.zeros((x.nrow(), y.nrow()), dtype=bool)
        if len(out[0]) > 0:
            rows = np.array(out[0], dtype=int)
            cols = np.array(out[1], dtype=int)
            m[rows, cols] = True
        return m


def relate_self(
    x: SpatVector,
    relation: str = "intersects",
    *,
    symmetrical: bool = True,
    na_rm: bool = True,
) -> np.ndarray:
    """
    Compute a spatial relation matrix within *x* (against itself).

    Parameters
    ----------
    x : SpatVector
    relation : str
    symmetrical : bool
        If True, return only the upper triangle (excluding diagonal).
    na_rm : bool

    Returns
    -------
    numpy.ndarray, shape (n, 2) of [id1, id2] (1-based) pairs.
    """
    n = x.nrow()
    out = x.related_between(x, relation, na_rm)
    messages(x, "relate_self")
    if len(out[0]) == 0:
        return np.empty((0, 2), dtype=int)
    rows = np.array(out[0], dtype=int) + 1
    cols = np.array(out[1], dtype=int) + 1
    pairs = np.column_stack([rows, cols])
    # Remove self-relation
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    if symmetrical:
        pairs = pairs[pairs[:, 0] < pairs[:, 1]]
    return pairs
