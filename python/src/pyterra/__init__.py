"""
pyterra — Python bindings for the terra spatial C++ library.

All classes are exposed directly from the compiled C++ extension
``_pyterra``.  The parallel object model mirrors the R terra package:

    SpatRaster            — raster data (one or more layers)
    SpatVector            — vector data (points, lines, polygons)
    SpatExtent            — a bounding box
    SpatOptions           — processing options (file format, memory, ...)
    SpatDataFrame         — attribute table
    SpatFactor            — categorical column
    SpatTime_v            — time dimension metadata
    SpatSRS               — spatial reference system helper
    SpatMessages          — error/warning message container
    SpatCategories        — raster category table
    SpatVectorCollection  — list of SpatVector objects
    SpatVectorProxy       — lazy-reading proxy for a SpatVector
    SpatRasterCollection  — list of SpatRaster objects
    SpatRasterStack       — multi-dataset raster stack (e.g. netCDF)
"""

from ._pyterra import (  # noqa: F401
    SpatRaster,
    SpatVector,
    SpatExtent,
    SpatOptions,
    SpatDataFrame,
    SpatFactor,
    SpatTime_v,
    SpatSRS,
    SpatMessages,
    SpatCategories,
    SpatVectorCollection,
    SpatVectorProxy,
    SpatRasterCollection,
    SpatRasterStack,
)

__all__ = [
    "SpatRaster",
    "SpatVector",
    "SpatExtent",
    "SpatOptions",
    "SpatDataFrame",
    "SpatFactor",
    "SpatTime_v",
    "SpatSRS",
    "SpatMessages",
    "SpatCategories",
    "SpatVectorCollection",
    "SpatVectorProxy",
    "SpatRasterCollection",
    "SpatRasterStack",
]
