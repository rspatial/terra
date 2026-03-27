"""
aggregate.py — spatial aggregation and disaggregation of raster and vector objects.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Union
import numpy as np

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# aggregate SpatRaster
# ---------------------------------------------------------------------------

_AGG_FUNS = {
    "sum", "mean", "median", "modal", "min", "max", "prod", "any", "all",
    "count", "sd", "std", "first",
}


def aggregate(
    x: SpatRaster,
    fact: Union[int, List[int]],
    fun: Union[str, Callable] = "mean",
    *,
    na_rm: bool = True,
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Aggregate cells of a SpatRaster.

    Parameters
    ----------
    x : SpatRaster
    fact : int or list of int
        Aggregation factor.  A scalar applies to rows, columns, and layers
        (if multi-layer).  A 2-element list [row_factor, col_factor]
        aggregates rows and columns separately.  A 3-element list
        [row_factor, col_factor, layer_factor] also aggregates layers.
    fun : str or callable
        Aggregation function.  Built-in shortcuts: ``"sum"``, ``"mean"``,
        ``"median"``, ``"min"``, ``"max"``, ``"sd"``, ``"modal"``, …
    na_rm : bool
        Ignore NA values.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    if isinstance(fact, int):
        fact_r, fact_c, fact_l = fact, fact, 1
    elif len(fact) == 2:
        fact_r, fact_c, fact_l = int(fact[0]), int(fact[1]), 1
    else:
        fact_r, fact_c, fact_l = int(fact[0]), int(fact[1]), int(fact[2])

    txt = fun if isinstance(fun, str) else getattr(fun, "__name__", "")
    if txt in _AGG_FUNS:
        opt = SpatOptions(filename, overwrite)
        xc = x.aggregate([fact_r, fact_c, fact_l], txt, na_rm, opt)
        return messages(xc, "aggregate")

    if not callable(fun):
        raise ValueError(f"Unknown aggregation function: {fun!r}")

    # Fallback: manual aggregation
    nr, nc = x.nrow(), x.ncol()
    nl = x.nlyr()
    vals = np.array(x.readValues(0, nr, 0, nc), dtype=float).reshape(nr * nc, nl, order='F')
    new_nr = max(1, nr // fact_r)
    new_nc = max(1, nc // fact_c)

    result = np.full((new_nr * new_nc, nl), float("nan"))
    for r in range(new_nr):
        for c in range(new_nc):
            r0, r1 = r * fact_r, min((r + 1) * fact_r, nr)
            c0, c1 = c * fact_c, min((c + 1) * fact_c, nc)
            cell_idx = []
            for ri in range(r0, r1):
                for ci in range(c0, c1):
                    cell_idx.append(ri * nc + ci)
            block = vals[cell_idx]
            if na_rm:
                block = block[~np.isnan(block).any(axis=1)] if block.ndim == 2 else block[~np.isnan(block)]
            if len(block) == 0:
                result[r * new_nc + c] = float("nan")
            else:
                result[r * new_nc + c] = fun(block, axis=0)

    from .rast import rast
    out = rast(x)
    out.setRows(new_nr)
    out.setCols(new_nc)
    opt = _opt()
    out.setValues(result.ravel(order='C').tolist(), opt)
    return out


def disagg(
    x: SpatRaster,
    fact: Union[int, List[int]],
    method: str = "near",
    filename: str = "",
    overwrite: bool = False,
) -> SpatRaster:
    """
    Disaggregate (increase resolution) of a SpatRaster.

    Parameters
    ----------
    x : SpatRaster
    fact : int or list of int
        Disaggregation factor.
    method : str
        ``"near"`` (nearest-neighbour, default) or ``"bilinear"``.
    filename : str
    overwrite : bool

    Returns
    -------
    SpatRaster
    """
    if isinstance(fact, int):
        fact = [fact, fact]
    elif len(fact) == 1:
        fact = [int(fact[0]), int(fact[0])]
    else:
        fact = [int(fact[0]), int(fact[1])]
    opt = SpatOptions(filename, overwrite)
    xc = x.disaggregate(fact, method, opt)
    return messages(xc, "disagg")


# ---------------------------------------------------------------------------
# aggregate SpatVector
# ---------------------------------------------------------------------------

def aggregate_vect(
    x: SpatVector,
    by: Optional[Union[str, List[str]]] = None,
    dissolve: bool = True,
) -> SpatVector:
    """
    Aggregate (dissolve) features of a SpatVector.

    Parameters
    ----------
    x : SpatVector
    by : str or list of str, optional
        Attribute column(s) to group by.
    dissolve : bool
        If True, merge geometries for each group.

    Returns
    -------
    SpatVector
    """
    if by is None:
        by = []
    elif isinstance(by, str):
        by = [by]
    opt = _opt()
    xc = x.aggregate(by, dissolve, opt)
    return messages(xc, "aggregate_vect")
