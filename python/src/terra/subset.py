"""
subset.py — layer and feature subsetting for SpatRaster and SpatVector.

``SpatRaster`` indexing (patched onto the pybind11 type) uses **0-based**
indices, like NumPy and typical Python APIs (unlike R's terra):

- **``r[[i]]``**, **``r[[i, j, ...]]``**, **``r[["name"]]``** — double-bracket
  style with a **list** key selects layers by **0-based** index or by name.
- **``r[i]``** — linear **cell** index (**0-based**); returns a 1-D ``numpy.ndarray``
  of values (one per layer).
- **``r[row, col]``** — **cell** at **0-based** row and column.

The low-level ``SpatRaster.subset(layers, opt)`` call still expects **0-based**
layer indices in the sequence (used internally by ``generics`` and others).
"""
from __future__ import annotations
import math
from typing import List, Optional, Sequence, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages, spatoptions


def _opt() -> SpatOptions:
    return SpatOptions()


def _is_scalar_int_index(x) -> bool:
    """True for a single integer index (not bool)."""
    if isinstance(x, bool):
        return False
    if isinstance(x, (int, np.integer)):
        return True
    return False


def _normalize_indices_0based(
    subset: Sequence,
    n: int,
    *,
    what: str = "index",
) -> List[int]:
    """
    Convert int / list of int to 0-based indices in ``[0, n)``.
    Negative ints count from the end (``-1`` is last). Bool list selects by position.
    """
    if len(subset) == 0:
        return []

    if isinstance(subset[0], (bool, np.bool_)):
        if len(subset) != n:
            raise ValueError(f"boolean mask length must be {n}, got {len(subset)}")
        return [idx for idx, v in enumerate(subset) if v]

    out: List[int] = []
    for v in subset:
        if v is None:
            continue
        j = int(v)
        if j < 0:
            j = n + j
        if j < 0 or j >= n:
            raise IndexError(f"{what} {v!r} out of range for length {n}")
        out.append(j)
    return out


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

    Equivalent to **``x[[subset]]``** (double brackets with a list key).
    Layer indices are **0-based** (``0 .. nlyr-1``); ``-1`` is the last layer.

    Parameters
    ----------
    x : SpatRaster
    subset : int, list of int, str, or list of str
        Layer indices (0-based) or layer names to select.
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
        indices: List[int] = []
        for name in subset:
            if name not in nms:
                raise ValueError(f"invalid layer name {name!r}")
            indices.append(nms.index(name))
        subset = indices

    if not isinstance(subset, (list, tuple, np.ndarray)):
        subset = [subset]

    pos = _normalize_indices_0based(list(subset), nl, what="layer")
    if negate:
        pos = [i for i in range(nl) if i not in set(pos)]

    if not pos:
        raise ValueError("subset results in no layers")

    opt = spatoptions(filename, overwrite)
    xc = x.subset(pos, opt)
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

    Row indices are **0-based** (``0 .. nrow-1``); ``-1`` is the last row.

    Parameters
    ----------
    x : SpatVector
    subset : list of bool or int, optional
        Row indices (0-based) or boolean mask selecting features.
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
        if hasattr(subset, "tolist"):
            subset = subset.tolist()
        if len(subset) > 0 and isinstance(subset[0], (bool, np.bool_)):
            row_idx = [i for i, v in enumerate(subset) if v]
        else:
            pos = _normalize_indices_0based(list(subset), x.nrow(), what="row")
            row_idx = pos

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


def _patch_spatraster_subset() -> None:
    """
    ``SpatRaster.subset(layers, opt)`` from pybind11 requires a sequence and
    ``SpatOptions``. Accept a single int (**0-based** layer index, same rules as
    :func:`subset_rast`) so ``r.subset(0)`` returns the first layer. Sequences
    are passed through as **0-based** layer indices for internal callers.
    """
    _cpp_subset = SpatRaster.subset

    def subset(self, layers, opt=None):
        if opt is None:
            opt = SpatOptions()
        if _is_scalar_int_index(layers):
            pos = _normalize_indices_0based([int(layers)], self.nlyr(), what="layer")
            if not pos:
                raise ValueError("subset results in no layers")
            return _cpp_subset(self, pos, opt)
        return _cpp_subset(self, layers, opt)

    SpatRaster.subset = subset  # type: ignore[assignment]


def _extract_cell_values(r: SpatRaster, cell_0based: float) -> np.ndarray:
    """Values at one cell (0-based index), one float per layer."""
    opt = SpatOptions()
    layers = r.extractCell([float(cell_0based)], opt)
    return np.array([row[0] for row in layers], dtype=float)


def _spatraster_getitem(self, key):
    """
    ``r[[layers...]]`` — *key* is ``list`` or 1-D array → layers (see :func:`subset_rast`).

    ``r[cell]`` — *key* is ``int`` → linear cell (**0-based**). ``r[row, col]`` — two ints.

    ``r[(i,)]`` is treated like ``r[i]``.
    """
    # --- Layer selection: r[[0]], r[[0,2]], r[["a"]]; key is list ---
    if isinstance(key, list):
        return subset_rast(self, key)

    # 1-D array: layer indices or mask
    if isinstance(key, np.ndarray):
        if key.ndim != 1:
            raise TypeError(
                "SpatRaster[[...]] expects a 1-D list or array of layer indices or names"
            )
        return subset_rast(self, key.tolist())

    if isinstance(key, slice) or key is Ellipsis:
        raise TypeError(
            "SpatRaster does not support slicing here; use r[[i]] for a layer, "
            "or r[i] / r[i, j] for cell values"
        )

    # --- (row, col) or (i,) ---
    if isinstance(key, tuple):
        if len(key) == 1:
            return _spatraster_getitem(self, key[0])
        if len(key) == 2 and all(_is_scalar_int_index(x) for x in key):
            r, c = int(key[0]), int(key[1])
            nr, nc = self.nrow(), self.ncol()
            if r < 0:
                r = nr + r
            if c < 0:
                c = nc + c
            if r < 0 or r >= nr or c < 0 or c >= nc:
                raise IndexError("row or column out of range")
            cell0 = self.cellFromRowCol([r], [c])[0]
            if math.isnan(cell0):
                raise IndexError("row or column out of range")
            return _extract_cell_values(self, cell0)
        raise TypeError(
            "SpatRaster indexing: use r[[layer]] for layers; "
            "r[cell] or r[row, col] for cell values (0-based)"
        )

    # --- Linear cell index (0-based) ---
    if _is_scalar_int_index(key):
        c = int(key)
        n = self.ncell()
        if c < 0:
            c = n + c
        if c < 0 or c >= n:
            raise IndexError("cell index out of range")
        return _extract_cell_values(self, float(c))

    raise TypeError(
        f"SpatRaster cannot index with {type(key).__name__!r}; "
        "use r[[...]] for layers, r[i] or r[i, j] for cells"
    )


def _patch_spatraster_getitem() -> None:
    SpatRaster.__getitem__ = _spatraster_getitem  # type: ignore[assignment]


_patch_spatraster_subset()
_patch_spatraster_getitem()
