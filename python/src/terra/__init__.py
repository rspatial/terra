"""
terra — Python bindings for the terra geospatial C++ library.

The compiled extension ``_terra`` exposes the same core C++ classes as the R
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
  plot()              pt.plot()
  plot_rgb()          pt.plot_rgb()
  ... and more in generics.py / arith.py / geom.py

Arithmetic operators (Arith_generics.R):
  r + r2, r - n, r * r2, etc.   operator overloads on SpatRaster
  e + e2 (union), e * n (scale)  operator overloads on SpatExtent
  v + v2 (union), v - v2 (erase) operator overloads on SpatVector
  is_na(), not_na(), which_max(), rast_sum(), compare_rast(), …

Geometry operations (geom.R):
  is_valid(), make_valid(), union_vect(), intersect_vect(), erase(),
  symdif(), buffer_vect(), hull(), voronoi(), delaunay(), spin(), …
"""

from ._terra import (  # noqa: F401
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
from .show import repr_extent, repr_raster, repr_vector, show, register_reprs  # noqa: F401
register_reprs()  # attach __repr__ / __str__ to C++ types
from .extent import ext                                               # noqa: F401
from .rast import rast                                                # noqa: F401
from .vect import vect                                                # noqa: F401
from .plot import plot, plot_rgb                                       # noqa: F401

from .arith import (                                                  # noqa: F401
    # NA / logical tests
    is_na, not_na, is_true, is_false,
    is_nan, is_finite, is_infinite,
    any_na, all_na, no_na, count_na,
    # summaries
    which_max, which_min, which_lyr,
    where_max, where_min,
    rast_sum, rast_mean, rast_min, rast_max,
    rast_median, rast_modal, stdev_rast,
    # compare / logic
    compare_rast, logic_rast_fn,
    # type coercion / queries
    as_int_rast, as_bool_rast,
    is_bool_rast, is_int_rast, is_num_rast,
    register_operators,
)
register_operators()  # attach +, -, *, /, ==, … to C++ types

from .geom import (                                                   # noqa: F401
    # validity
    is_valid, make_valid,
    # set operations
    union_vect, intersect_vect, erase, symdif, cover_vect,
    # crop / mask
    crop_vect, mask_vect,
    # geometry modifications
    buffer_vect, disagg_vect, flip_vect, spin,
    hull, delaunay, voronoi, elongate,
    merge_lines, make_nodes, remove_dup_nodes,
    simplify_geom, thin_geom,
    shared_paths, snap_vect, gaps,
    force_ccw, width_vect, clearance,
    # predicates
    is_empty,
)

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
    "plot", "plot_rgb",
    "messages", "character_crs",
    "show", "repr_raster", "repr_vector", "repr_extent",
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
    # arith (Arith_generics.R)
    "is_na", "not_na", "is_true", "is_false",
    "is_nan", "is_finite", "is_infinite",
    "any_na", "all_na", "no_na", "count_na",
    "which_max", "which_min", "which_lyr",
    "where_max", "where_min",
    "rast_sum", "rast_mean", "rast_min", "rast_max",
    "rast_median", "rast_modal", "stdev_rast",
    "compare_rast", "logic_rast_fn",
    "as_int_rast", "as_bool_rast",
    "is_bool_rast", "is_int_rast", "is_num_rast",
    # geom (geom.R)
    "is_valid", "make_valid",
    "union_vect", "intersect_vect", "erase", "symdif", "cover_vect",
    "crop_vect", "mask_vect",
    "buffer_vect", "disagg_vect", "flip_vect", "spin",
    "hull", "delaunay", "voronoi", "elongate",
    "merge_lines", "make_nodes", "remove_dup_nodes",
    "simplify_geom", "thin_geom",
    "shared_paths", "snap_vect", "gaps",
    "force_ccw", "width_vect", "clearance",
    "is_empty",
    # Core types (C++)
    "SpatRaster", "SpatVector", "SpatExtent", "SpatOptions",
    "SpatDataFrame", "SpatFactor", "SpatTime_v", "SpatSRS",
    "SpatMessages", "SpatCategories", "SpatVectorCollection",
    "SpatVectorProxy", "SpatRasterCollection", "SpatRasterStack",
    "__version__",
]
