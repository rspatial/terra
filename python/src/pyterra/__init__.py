"""
pyterra — Python bindings for the terra spatial C++ library.

The compiled extension ``_pyterra`` exposes the same core C++ classes as the R
package.  The functions in this package mirror the **R** ``terra`` API so that
workflows can be translated with minimal renaming.

Quick reference (R → Python):
  rast()              pt.rast()
  vect()              pt.vect()
  ext()               pt.ext()
  crs()               pt.crs()
  nrow/ncol/nlyr      pt.nrow / pt.ncol / pt.nlyr
  res()               pt.res()
  origin()            pt.origin()
  crop()              pt.crop()
  mask()              pt.mask()
  project(rast)       pt.project_raster()
  project(vect)       pt.project_vector()
  resample()          pt.resample()
  terrain()           pt.terrain()
  classify()          pt.classify()
  flip()              pt.flip()
  trim()              pt.trim()
  boundaries()        pt.boundaries()
  patches()           pt.patches()
  cellSize()          pt.cellSize()
  ... and more in generics.py
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

from ._helpers import character_crs, messages                        # noqa: F401
from .crs import crs                                                  # noqa: F401
from .extent import ext                                               # noqa: F401
from .rast import rast                                                # noqa: F401
from .vect import vect                                                # noqa: F401
from .generics import (                                               # noqa: F401
    # dimensions / metadata
    nrow, ncol, nlyr, ncell, res, origin,
    # helpers
    spat_options, deepcopy, tighten,
    # extent
    ext_align,
    # raster geometry
    is_rotated, is_flipped, flip, rotate, shift, rescale,
    trans, trim, rev_raster,
    # raster values
    clamp, clamp_ts, classify, subst, cover, diff_raster,
    disagg, segregate, selectRange, sort_raster,
    range_fill, weighted_mean,
    # raster analysis
    boundaries, patches, cellSize, surfArea, terrain,
    sieve, rectify, stretch, scale_linear, scale_raster,
    quantile_raster, atan_2,
    # raster processing
    crop, mask, project_raster, resample, intersect_rast,
    # vector
    project_vector, shift_vect, rotate_vect, rescale_vect, trans_vect,
    # scoff
    scoff, scoff_set,
)

__version__ = "0.1.0"

__all__ = [
    # High-level API (R-like)
    "rast", "vect", "ext", "crs",
    "messages", "character_crs",
    # dimensions
    "nrow", "ncol", "nlyr", "ncell", "res", "origin",
    # helpers
    "spat_options", "deepcopy", "tighten",
    # extent
    "ext_align",
    # raster geometry
    "is_rotated", "is_flipped", "flip", "rotate", "shift", "rescale",
    "trans", "trim", "rev_raster",
    # raster values
    "clamp", "clamp_ts", "classify", "subst", "cover", "diff_raster",
    "disagg", "segregate", "selectRange", "sort_raster",
    "range_fill", "weighted_mean",
    # raster analysis
    "boundaries", "patches", "cellSize", "surfArea", "terrain",
    "sieve", "rectify", "stretch", "scale_linear", "scale_raster",
    "quantile_raster", "atan_2",
    # raster processing
    "crop", "mask", "project_raster", "resample", "intersect_rast",
    # vector
    "project_vector", "shift_vect", "rotate_vect", "rescale_vect", "trans_vect",
    # scoff
    "scoff", "scoff_set",
    # Core types (C++)
    "SpatRaster", "SpatVector", "SpatExtent", "SpatOptions",
    "SpatDataFrame", "SpatFactor", "SpatTime_v", "SpatSRS",
    "SpatMessages", "SpatCategories", "SpatVectorCollection",
    "SpatVectorProxy", "SpatRasterCollection", "SpatRasterStack",
    "__version__",
]
