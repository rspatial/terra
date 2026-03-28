"""
names.py — layer and variable names for SpatRaster and SpatVector.
"""
from __future__ import annotations
from typing import List, Optional, Union

from ._terra import SpatRaster, SpatVector, SpatOptions


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# SpatRaster layer names
# ---------------------------------------------------------------------------

def names_rast(x: SpatRaster) -> List[str]:
    """Return the layer names of *x*."""
    return list(x.names)


def set_names_rast(
    x: SpatRaster,
    value: List[str],
    index: Optional[List[int]] = None,
    validate: bool = False,
) -> SpatRaster:
    """
    Return a copy of *x* with new layer names.

    Parameters
    ----------
    x : SpatRaster
    value : list of str
        New names.  Length must equal nlyr(x) (or len(index) if provided).
    index : list of int, optional
        0-based layer indices to rename.  If None, rename all layers.
    validate : bool
        If True, sanitize names to be valid identifiers.

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    if index is None:
        index = list(range(xc.nlyr()))
    if len(value) != len(index):
        raise ValueError("length of value does not match length of index")
    current = list(xc.names)
    for i, v in zip(index, value):
        current[i] = str(v)
    if validate:
        import re
        seen: dict = {}
        valid = []
        for n in current:
            n = re.sub(r'[^A-Za-z0-9_.]', '.', n)
            if not n or n[0].isdigit():
                n = 'X' + n
            base = n
            cnt = seen.get(base, 0)
            if cnt:
                n = f"{base}.{cnt}"
            seen[base] = cnt + 1
            valid.append(n)
        current = valid
    if not xc.setNames(current, False):
        raise RuntimeError("cannot set these names")
    return xc


def set_names_inplace(
    x: SpatRaster,
    value: List[str],
    index: Optional[List[int]] = None,
    validate: bool = False,
) -> None:
    """
    Modify layer names of *x* in-place (no copy).

    Parameters
    ----------
    x : SpatRaster
    value : list of str
    index : list of int, optional
        0-based layer indices.
    validate : bool
    """
    if index is None:
        index = list(range(x.nlyr()))
    current = list(x.names)
    for i, v in zip(index, value):
        current[i] = str(v)
    if not x.setNames(current, False):
        raise RuntimeError("cannot set these names")


# ---------------------------------------------------------------------------
# SpatVector attribute names
# ---------------------------------------------------------------------------

def names_vect(x: SpatVector) -> List[str]:
    """Return the attribute column names of *x*."""
    return list(x.names)


def set_names_vect(x: SpatVector, value: List[str]) -> SpatVector:
    """
    Return a copy of *x* with new attribute column names.

    Parameters
    ----------
    x : SpatVector
    value : list of str
        New column names.  Length must match ncol(x).

    Returns
    -------
    SpatVector
    """
    if len(value) != x.ncol():
        raise ValueError("incorrect number of names")
    xc = x.deepcopy()
    xc.names = [str(v) for v in value]
    return xc


# ---------------------------------------------------------------------------
# varnames — source/variable names embedded in a file
# ---------------------------------------------------------------------------

def varnames(x: SpatRaster) -> List[str]:
    """Return the variable (source) names of *x*."""
    return list(x.get_sourcenames())


def set_varnames(x: SpatRaster, value: List[str]) -> SpatRaster:
    """
    Return a copy of *x* with new variable names.

    Parameters
    ----------
    x : SpatRaster
    value : list of str

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    if not xc.set_sourcenames([str(v) for v in value]):
        raise RuntimeError("cannot set varnames")
    return xc


# ---------------------------------------------------------------------------
# longnames — human-readable long layer names
# ---------------------------------------------------------------------------

def longnames(x: SpatRaster) -> List[str]:
    """Return the long names of *x*."""
    return list(x.get_sourcenames_long())


def set_longnames(x: SpatRaster, value: List[str]) -> SpatRaster:
    """
    Return a copy of *x* with new long names.

    Parameters
    ----------
    x : SpatRaster
    value : list of str

    Returns
    -------
    SpatRaster
    """
    xc = x.deepcopy() if hasattr(x, 'deepcopy') else x
    if not xc.set_sourcenames_long([str(v) for v in value]):
        raise RuntimeError("cannot set longnames")
    return xc
