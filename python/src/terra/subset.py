"""
subset.py — layer and feature subsetting for SpatRaster and SpatVector.
"""
from __future__ import annotations
from typing import List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


def _positive_indices(i, n: int) -> List[int]:
    """Convert an index list (1-based, can include negatives) to positive 1-based ints."""
    if isinstance(i, (bool, np.bool_)):
        i = [i]
    if hasattr(i, 'tolist'):
        i = i.tolist()
    if not isinstance(i, (list, tuple)):
        i = [i]

    if len(i) == 0:
        return []

    if isinstance(i[0], (bool, np.bool_)):
        return [idx + 1 for idx, v in enumerate(i) if v]

    all_neg = all(v <= 0 for v in i if v is not None)
    all_pos = all(v >= 0 for v in i if v is not None)
    if not (all_neg or all_pos):
        raise ValueError("cannot mix positive and negative indices")

    all_idx = list(range(1, n + 1))
    if all_neg:
        exclude = {abs(v) for v in i if v is not None}
        return [idx for idx in all_idx if idx not in exclude]
    return [int(v) for v in i if v is not None and 1 <= int(v) <= n]


# ---------------------------------------------------------------------------
# SpatRaster subset
# ---------------------------------------------------------------------------

def subset_rast(
    x: SpatRaster,
    subset: Union[int, List, str, List[str]],
    *,
    negate: bool = False,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Select a subset of layers from *x*.

    Parameters
    ----------
    x : SpatRaster
    subset : int, list of int, str, or list of str
        Layer indices (1-based) or layer names to select.
    negate : bool
        If True, drop the selected layers instead.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    nl = x.nlyr()
    nms = list(x.names)

    if isinstance(subset, str):
        subset = [subset]
    if isinstance(subset, (list, tuple)) and len(subset) > 0 and isinstance(subset[0], str):
        indices = []
        for name in subset:
            if name not in nms:
                raise ValueError(f"invalid layer name {name!r}")
            indices.append(nms.index(name) + 1)
        subset = indices

    if not isinstance(subset, (list, tuple, np.ndarray)):
        subset = [subset]

    pos = _positive_indices(subset, nl)
    if negate:
        pos = [i for i in range(1, nl + 1) if i not in pos]

    if not pos:
        raise ValueError("subset results in no layers")

    opt = SpatOptions(filename, overwrite)
    xc = x.subset([i - 1 for i in pos], opt)
    return messages(xc, "subset")


# ---------------------------------------------------------------------------
# SpatVector subset (by rows and columns)
# ---------------------------------------------------------------------------

def subset_vect(
    x: SpatVector,
    subset: Optional[Union[List[bool], List[int], np.ndarray]] = None,
    select: Optional[Union[str, List[str]]] = None,
) -> SpatVector:
    """
    Select features and/or columns from a SpatVector.

    Parameters
    ----------
    x : SpatVector
    subset : list of bool or int, optional
        Row indices (1-based) or boolean mask selecting features.
    select : str or list of str, optional
        Column names to keep.

    Returns
    -------
    SpatVector
    """
    # Row subsetting
    if subset is None:
        row_idx = list(range(x.nrow()))
    else:
        if hasattr(subset, 'tolist'):
            subset = subset.tolist()
        if len(subset) > 0 and isinstance(subset[0], (bool, np.bool_)):
            row_idx = [i for i, v in enumerate(subset) if v]
        else:
            pos = _positive_indices(subset, x.nrow())
            row_idx = [i - 1 for i in pos]

    xc = x.subset(row_idx, _opt())

    # Column subsetting
    if select is not None:
        if isinstance(select, str):
            select = [select]
        all_cols = list(xc.names)
        keep_idx = []
        for col in select:
            if col not in all_cols:
                raise ValueError(f"column {col!r} not found")
            keep_idx.append(all_cols.index(col))
        xc = xc.selectRange(keep_idx)

    return messages(xc, "subset_vect")
