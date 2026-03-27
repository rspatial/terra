"""
terra — Python bindings for the terra geospatial C++ library.

The compiled extension ``_terra`` exposes the same core C++ classes as the R
package.  The functions in this package mirror the **R** ``terra`` API so that
workflows can be translated with minimal renaming.

Two equivalent calling styles are supported:

  Functional style (R-like):       Method style (Pythonic):
  ─────────────────────────────    ─────────────────────────────────
  pt.crop(r, e)                    r.crop(e)
  pt.mask(r, mask_r)               r.mask(mask_r)
  pt.aggregate(r, 2)               r.aggregate(2)
  pt.focal(r, 3)                   r.focal(3)
  pt.values(r)                     r.values()
  pt.buffer_vect(v, 1000)          v.buffer(1000)
  pt.project_vector(v, crs)        v.project(crs)

Quick reference (R → Python):
  rast()              pt.rast()
  vect()              pt.vect()
  ext()               pt.ext()
  crs()               pt.crs()
  crop()              pt.crop(r, e)  or  r.crop(e)
  mask()              pt.mask(r, m)  or  r.mask(m)
  project(rast)       pt.project_raster()  or  r.project(crs)
  project(vect)       pt.project_vector()  or  v.project(crs)
  resample()          pt.resample()  or  r.resample(template)
  classify()          pt.classify()  or  r.classify(rcl)
  terrain()           pt.terrain()  or  r.terrain()
  focal()             pt.focal()  or  r.focal(w)
  aggregate()         pt.aggregate()  or  r.aggregate(fact)
  values()            pt.values(r)  or  r.values()
  plot()              pt.plot(r)  or  r.plot()
  ... and more in generics.py / arith.py / geom.py / methods.py

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

from .methods import register_methods                                 # noqa: F401
register_methods()    # attach r.crop(), r.mask(), v.buffer(), … to C++ types

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

# ---- New translation modules -----------------------------------------------
from .values import (                                                 # noqa: F401
    has_values, in_memory, sources,
    has_min_max, min_max, set_min_max,
    values, set_values, focal_values,
    vect_values, set_vect_values,
    compare_geom,
)
from .levels import (                                                 # noqa: F401
    is_factor, as_factor,
    levels, set_levels,
    cats, set_cats, categories,
    active_cat, set_active_cat,
    add_cats, drop_levels, concats, catalyze,
    has_colors, coltab, set_coltab,
)
from .names import (                                                  # noqa: F401
    names_rast, set_names_rast, set_names_inplace,
    names_vect, set_names_vect,
    varnames, set_varnames,
    longnames, set_longnames,
)
from .app import app, lapp, tapp, xapp, rapp, sapp                   # noqa: F401
from .focal import focal, focal3D, focal_mat                          # noqa: F401
from .aggregate import aggregate, disagg as aggregate_disagg, aggregate_vect  # noqa: F401
from .zonal import zonal                                              # noqa: F401
from .crosstab import crosstab                                        # noqa: F401
from .freq import freq                                                # noqa: F401
from .flow_accumulation import flow_accumulation                     # noqa: F401
from .pitfinder import pitfinder                                     # noqa: F401
from .extract import extract, extract_xy                              # noqa: F401
from .math import (                                                   # noqa: F401
    math, log, sqrt, abs_ as rast_abs, ceiling, floor,
    round_, cumsum, cumprod, cummax, cummin,
    floor_ext, ceiling_ext, round_ext,
    ifel,
)
from .cells import (                                                  # noqa: F401
    cells,
    row_from_y, col_from_x,
    cell_from_xy, cell_from_row_col,
    xy_from_cell, row_col_from_cell,
)
from .init import init                                                # noqa: F401
from .distance import (                                               # noqa: F401
    buffer_rast, distance_rast,
    cost_dist, grid_dist,
    distance_xy, distance_vect_self, distance_vect, distance_points,
)
from .rasterize import rasterize, rasterize_geom                     # noqa: F401
from .time import has_time, time_info, get_time, set_time            # noqa: F401
from .write import (                                                  # noqa: F401
    write_raster, write_start, write_values, write_stop, blocks,
    write_vector,
)
from .sample import spat_sample, grid_sample                         # noqa: F401
from .stats import (                                                  # noqa: F401
    row_sums, col_sums, row_means, col_means,
    match_rast, is_in,
    autocor, layer_cor,
)
from .merge import merge as merge_rast, mosaic, merge_vect           # noqa: F401
from .relate import is_related, relate, relate_self                  # noqa: F401
from .subset import subset_rast, subset_vect                         # noqa: F401
from .window import has_window, set_window, remove_window, extend    # noqa: F401
from .coerce import (                                                 # noqa: F401
    as_polygons, as_lines, as_points,
    as_array, as_matrix, as_data_frame,
)
from .spatvec import (                                                # noqa: F401
    geomtype, is_lines, is_polygons, is_points,
    geom, crds,
    expanse, perim, nseg,
    fill_holes, vect_as_df, geom_as_wkt,
)

__version__ = "0.1.0"

__all__ = [
    # High-level API (R-like)
    "rast", "vect", "ext", "crs",
    "register_methods",
    "plot", "plot_rgb",
    "messages", "character_crs",
    "show", "repr_raster", "repr_vector", "repr_extent",
    # values
    "has_values", "in_memory", "sources",
    "has_min_max", "min_max", "set_min_max",
    "values", "set_values", "focal_values",
    "vect_values", "set_vect_values", "compare_geom",
    # levels / colors
    "is_factor", "as_factor",
    "levels", "set_levels",
    "cats", "set_cats", "categories",
    "active_cat", "set_active_cat",
    "add_cats", "drop_levels", "concats", "catalyze",
    "has_colors", "coltab", "set_coltab",
    # names
    "names_rast", "set_names_rast", "set_names_inplace",
    "names_vect", "set_names_vect",
    "varnames", "set_varnames",
    "longnames", "set_longnames",
    # app
    "app", "lapp", "tapp", "xapp", "rapp", "sapp",
    # focal
    "focal", "focal3D", "focal_mat",
    # aggregate
    "aggregate", "aggregate_disagg", "aggregate_vect",
    # zonal
    "zonal",
    # crosstab
    "crosstab",
    # freq
    "freq",
    # flow accumulation
    "flow_accumulation",
    # pitfinder
    "pitfinder",
    # extract
    "extract", "extract_xy",
    # math
    "math", "log", "sqrt", "rast_abs", "ceiling", "floor",
    "round_", "cumsum", "cumprod", "cummax", "cummin",
    "floor_ext", "ceiling_ext", "round_ext",
    "ifel",
    # cells
    "cells",
    "row_from_y", "col_from_x",
    "cell_from_xy", "cell_from_row_col",
    "xy_from_cell", "row_col_from_cell",
    # init
    "init",
    # distance
    "buffer_rast", "distance_rast",
    "cost_dist", "grid_dist",
    "distance_xy", "distance_vect_self", "distance_vect", "distance_points",
    # rasterize
    "rasterize", "rasterize_geom",
    # time
    "has_time", "time_info", "get_time", "set_time",
    # write
    "write_raster", "write_start", "write_values", "write_stop", "blocks",
    "write_vector",
    # sample
    "spat_sample", "grid_sample",
    # stats
    "row_sums", "col_sums", "row_means", "col_means",
    "match_rast", "is_in",
    "autocor", "layer_cor",
    # merge
    "merge_rast", "mosaic", "merge_vect",
    # relate
    "is_related", "relate", "relate_self",
    # subset
    "subset_rast", "subset_vect",
    # window
    "has_window", "set_window", "remove_window", "extend",
    # coerce
    "as_polygons", "as_lines", "as_points",
    "as_array", "as_matrix", "as_data_frame",
    # spatvec
    "geomtype", "is_lines", "is_polygons", "is_points",
    "geom", "crds",
    "expanse", "perim", "nseg",
    "fill_holes", "vect_as_df", "geom_as_wkt",
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
