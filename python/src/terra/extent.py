"""Extent helpers — parallel to R ``terra::ext()``."""

from __future__ import annotations

from typing import Any

from ._helpers import messages
from ._terra import SpatExtent

__all__ = ["ext"]


def ext(
    *args: Any,
    xy: bool = False,
) -> SpatExtent:
    """
    Create or copy a :class:`SpatExtent` (bounding box), like R ``terra::ext()``.

    * ``ext()`` — empty extent object.
    * ``ext(xmin, xmax, ymin, ymax)`` — four numbers.
    * ``ext(SpatExtent)`` — deep copy.
    * ``ext([xmin, xmax, ymin, ymax])`` — length-4 sequence.
    * With **numpy**, a 2-column coordinate matrix uses column-wise min/max.
    """
    if len(args) == 0:
        e = SpatExtent()
        return messages(e, "ext")

    if len(args) == 1:
        a = args[0]
        if isinstance(a, SpatExtent):
            return messages(a.deepcopy(), "ext")

        try:
            import numpy as np  # type: ignore

            if isinstance(a, np.ndarray) and a.ndim == 2 and a.shape[1] == 2:
                mn = a.min(axis=0)
                mx = a.max(axis=0)
                e = SpatExtent(float(mn[0]), float(mx[0]), float(mn[1]), float(mx[1]))
                if not e.valid:
                    raise ValueError("ext: invalid extent")
                return messages(e, "ext")
        except ImportError:
            pass

        if isinstance(a, (list, tuple)) or (
            hasattr(a, "__iter__") and not isinstance(a, (str, bytes))
        ):
            v = list(a)
            if len(v) != 4:
                raise ValueError("ext: expected four numbers")
            x0, x1, x2, x3 = float(v[0]), float(v[1]), float(v[2]), float(v[3])
            if xy:
                e = SpatExtent(x0, x2, x1, x3)  # xmin, ymin, xmax, ymax
            else:
                e = SpatExtent(x0, x1, x2, x3)  # xmin, xmax, ymin, ymax
            if not e.valid:
                raise ValueError("ext: invalid extent")
            return messages(e, "ext")

        raise TypeError("ext: expected four numbers, SpatExtent, or a 2-column array")

    if len(args) == 4:
        x0, x1, x2, x3 = (
            float(args[0]),
            float(args[1]),
            float(args[2]),
            float(args[3]),
        )
        if xy:
            e = SpatExtent(x0, x2, x1, x3)
        else:
            e = SpatExtent(x0, x1, x2, x3)
        if not e.valid:
            raise ValueError("ext: invalid extent")
        return messages(e, "ext")

    raise TypeError("ext: invalid arguments")
