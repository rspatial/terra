"""
Display helpers for core terra objects.

Printing and ``__repr__`` use the same formatted strings as R (``show.R``) and
the C++ ``show()`` implementations in ``src/show.cpp`` — no duplicated layout
logic in Python.
"""

from __future__ import annotations

from typing import Any, Type, Tuple

__all__ = [
    "repr_extent",
    "repr_raster",
    "repr_vector",
    "show",
    "register_reprs",
]


def _show_text(obj: Any) -> str:
    """Single multi-line string from C++ (trailing newline optional)."""
    return obj.show().rstrip("\n")


def repr_extent(e: Any) -> str:
    return _show_text(e)


def repr_raster(r: Any) -> str:
    return _show_text(r)


def repr_vector(v: Any) -> str:
    return _show_text(v)


def show(obj: Any) -> None:
    """Print the same summary as R ``show()`` / ``print()`` for supported types."""
    from ._terra import (
        SpatExtent,
        SpatRaster,
        SpatRasterCollection,
        SpatRasterStack,
        SpatVector,
        SpatVectorCollection,
        SpatVectorProxy,
    )

    _WITH_SHOW: Tuple[Type[Any], ...] = (
        SpatExtent,
        SpatRaster,
        SpatVector,
        SpatRasterStack,
        SpatRasterCollection,
        SpatVectorCollection,
        SpatVectorProxy,
    )
    if isinstance(obj, _WITH_SHOW):
        print(obj.show(), end="")
    else:
        print(repr(obj))


def register_reprs() -> None:
    """
    Attach ``__repr__`` and ``__str__`` to the pybind11 types so the REPL shows
    the C++ ``show()`` text automatically.

    Called once at import time from ``__init__.py``.
    """
    from ._terra import (
        SpatExtent,
        SpatRaster,
        SpatRasterCollection,
        SpatRasterStack,
        SpatVector,
        SpatVectorCollection,
        SpatVectorProxy,
    )

    def _repr_str(self: Any) -> str:
        return _show_text(self)

    for cls in (
        SpatExtent,
        SpatRaster,
        SpatVector,
        SpatRasterStack,
        SpatRasterCollection,
        SpatVectorCollection,
        SpatVectorProxy,
    ):
        cls.__repr__ = _repr_str  # type: ignore[assignment]
        cls.__str__ = _repr_str  # type: ignore[assignment]
