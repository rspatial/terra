"""
methods.py — register high-level Python functions as methods on the C++ types.

After :func:`register_methods` is called (done once at import time by
``terra/__init__.py``), the Python-level functions from ``generics.py``,
``values.py``, ``cells.py``, etc. are accessible directly on objects:

    r.crop(e)           # instead of terra.crop(r, e)
    r.mask(mask_r)      # instead of terra.mask(r, mask_r)
    r.values()          # instead of terra.values(r)
    r.aggregate(2)      # instead of terra.aggregate(r, 2)
    v.buffer(1000)      # instead of terra.buffer_vect(v, 1000)
    e.intersect(e2)     # SpatExtent.intersect is already C++; this is the Python alias

Only the **Python-level** wrappers are registered here.  The raw C++ methods
(e.g. ``r.classify(...)``) remain accessible but may require different argument
conventions — the registered wrappers follow the R-like Python API.
"""
from __future__ import annotations

from ._terra import SpatRaster, SpatVector, SpatExtent


def register_methods() -> None:
    """
    Attach high-level Python functions as methods on SpatRaster, SpatVector,
    and SpatExtent.

    Called once at import time by ``terra/__init__.py``.
    """

    # ── delayed imports to avoid circular deps ────────────────────────────────
    # We import inside the function so this module can be imported at any time
    # without triggering circular import issues.

    # ── SpatRaster methods ────────────────────────────────────────────────────

    from .generics import (
        crop, mask as mask_rast, resample, project_raster,
        classify, subst, clamp, clamp_ts, cover, diff_raster,
        boundaries, patches, cellSize, surfArea, terrain,
        sieve, stretch, scale_linear, scale_raster,
        trim, flip, rotate, shift, rescale, rev_raster, sort_raster,
        disagg, segregate, selectRange, range_fill, weighted_mean,
        quantile_raster, atan_2,
    )
    from .values import (
        values, set_values, focal_values,
        has_values, in_memory, sources,
        has_min_max, min_max, set_min_max, compare_geom,
    )
    from .cells import (
        cells, cell_from_xy, cell_from_row_col,
        xy_from_cell, row_col_from_cell,
        row_from_y, col_from_x,
    )
    from .aggregate import aggregate
    from .focal import focal, focal3D, focal_mat
    from .zonal import zonal
    from .extract import extract
    from .app import app, lapp, tapp, sapp
    from .math import (
        math as rast_math, log, sqrt,
        abs_ as rast_abs, ceiling, floor, round_,
        cumsum, cumprod, cummax, cummin,
        ifel,
    )
    from .distance import distance_rast, cost_dist, grid_dist
    from .rasterize import rasterize
    from .sample import spat_sample
    from .window import has_window, set_window, remove_window, extend
    from .init import init
    from .arith import (
        is_na, not_na, any_na, all_na, no_na, count_na,
        is_nan, is_finite, is_infinite,
        which_max, which_min, which_lyr,
        where_max, where_min,
        rast_sum, rast_mean, rast_min, rast_max,
        rast_median, rast_modal, stdev_rast,
        compare_rast,
    )
    from .levels import (
        is_factor, as_factor,
        levels, set_levels,
        cats, set_cats, drop_levels, catalyze,
        has_colors, coltab, set_coltab,
    )
    from .names import (
        names_rast, set_names_rast, set_names_inplace,
        varnames, set_varnames, longnames, set_longnames,
    )
    from .time import has_time, time_info, get_time, set_time
    from .write import write_raster
    from .stats import row_sums, col_sums, row_means, col_means, autocor, layer_cor
    from .merge import merge as merge_rast, mosaic
    from .coerce import as_polygons, as_lines, as_points, as_array, as_matrix, as_data_frame
    from .plot import plot, plot_rgb

    _rast_methods = {
        # extent (R-style; backed by SpatExtent.vector = xmin, xmax, ymin, ymax)
        "xmin":           lambda self: float(self.extent.vector[0]),
        "xmax":           lambda self: float(self.extent.vector[1]),
        "ymin":           lambda self: float(self.extent.vector[2]),
        "ymax":           lambda self: float(self.extent.vector[3]),
        # geometry / metadata
        "crop":           lambda self, y, **kw: crop(self, y, **kw),
        "mask":           lambda self, m, **kw: mask_rast(self, m, **kw),
        "resample":       lambda self, y, **kw: resample(self, y, **kw),
        "project":        lambda self, crs, **kw: project_raster(self, crs, **kw),
        "trim":           lambda self, **kw: trim(self, **kw),
        "flip":           lambda self, **kw: flip(self, **kw),
        "rotate":         lambda self, **kw: rotate(self, **kw),
        "shift":          lambda self, **kw: shift(self, **kw),
        "rev":            lambda self, **kw: rev_raster(self, **kw),
        "aggregate":      lambda self, fact, fun="mean", **kw: aggregate(self, fact, fun, **kw),
        "disagg":         lambda self, fact, **kw: disagg(self, fact, **kw),
        # values
        "values":         lambda self: values(self),
        "set_values":     lambda self, v: set_values(self, v),
        "focal_values":   lambda self, w=3, **kw: focal_values(self, w, **kw),
        "has_values":     lambda self: has_values(self),
        "in_memory":      lambda self: in_memory(self),
        "min_max":        lambda self: min_max(self),
        # classification / recoding
        "classify":       lambda self, rcl, **kw: classify(self, rcl, **kw),
        "subst":          lambda self, from_v, to_v, **kw: subst(self, from_v, to_v, **kw),
        "clamp":          lambda self, lower, upper, **kw: clamp(self, lower, upper, **kw),
        "ifel":           lambda self, test, false_val, **kw: ifel(self, test, false_val, **kw),
        # raster analysis
        "boundaries":     lambda self, **kw: boundaries(self, **kw),
        "patches":        lambda self, **kw: patches(self, **kw),
        "terrain":        lambda self, v="slope", **kw: terrain(self, v, **kw),
        "focal":          lambda self, w, fun="sum", **kw: focal(self, w, fun, **kw),
        "focal3D":        lambda self, w, fun="sum", **kw: focal3D(self, w, fun, **kw),
        "zonal":          lambda self, z, fun="mean", **kw: zonal(self, z, fun, **kw),
        "app":            lambda self, fun, **kw: app(self, fun, **kw),
        "lapp":           lambda self, fun, **kw: lapp(self, fun, **kw),
        "tapp":           lambda self, index, fun, **kw: tapp(self, index, fun, **kw),
        "sapp":           lambda self, fun, **kw: sapp(self, fun, **kw),
        "extract":        lambda self, y, **kw: extract(self, y, **kw),
        "rasterize":      lambda self, v, **kw: rasterize(v, self, **kw),
        # cells / coordinates
        "cells":          lambda self, y=None, **kw: cells(self, y, **kw),
        "xy_from_cell":   lambda self, cell: xy_from_cell(self, cell),
        "cell_from_xy":   lambda self, xy: cell_from_xy(self, xy),
        # math
        "log":            lambda self, base=None: log(self) if base is None else log(self, base),
        "sqrt":           lambda self: sqrt(self),
        "abs":            lambda self: rast_abs(self),
        "ceiling":        lambda self: ceiling(self),
        "floor":          lambda self: floor(self),
        "round":          lambda self, digits=0: round_(self, digits),
        "cumsum":         lambda self: cumsum(self),
        "stretch":        lambda self, **kw: stretch(self, **kw),
        # NA
        "is_na":          lambda self: is_na(self),
        "not_na":         lambda self: not_na(self),
        "count_na":       lambda self: count_na(self),
        # summaries (global)
        "sum":            lambda self, na_rm=False: rast_sum(self, na_rm=na_rm),
        "mean":           lambda self, na_rm=False: rast_mean(self, na_rm=na_rm),
        "min":            lambda self, na_rm=False: rast_min(self, na_rm=na_rm),
        "max":            lambda self, na_rm=False: rast_max(self, na_rm=na_rm),
        "sd":             lambda self, na_rm=False: stdev_rast(self, na_rm=na_rm),
        # levels / categories
        "levels":         lambda self: levels(self),
        "set_levels":     lambda self, v: set_levels(self, v),
        "cats":           lambda self: cats(self),
        "is_factor":      lambda self: is_factor(self),
        "as_factor":      lambda self: as_factor(self),
        "catalyze":       lambda self: catalyze(self),
        # names
        "names_rast":     lambda self: names_rast(self),
        # time
        "has_time":       lambda self: has_time(self),
        "get_time":       lambda self: get_time(self),
        "set_time":       lambda self, v, **kw: set_time(self, v, **kw),
        # write
        "write":          lambda self, filename, **kw: write_raster(self, filename, **kw),
        # stats
        "autocor":        lambda self, **kw: autocor(self, **kw),
        "layer_cor":      lambda self, fun="cor", **kw: layer_cor(self, fun, **kw),
        "row_sums":       lambda self, **kw: row_sums(self, **kw),
        "col_sums":       lambda self, **kw: col_sums(self, **kw),
        # merge
        "merge":          lambda self, *others, **kw: merge_rast(self, *others, **kw),
        "mosaic":         lambda self, *others, **kw: mosaic(self, *others, **kw),
        # sampling
        "sample":         lambda self, size, **kw: spat_sample(self, size, **kw),
        # window
        "window":         lambda self: has_window(self),
        "set_window":     lambda self, e: set_window(self, e),
        "remove_window":  lambda self: remove_window(self),
        "extend":         lambda self, y, **kw: extend(self, y, **kw),
        # coerce
        "as_polygons":    lambda self, **kw: as_polygons(self, **kw),
        "as_points":      lambda self, **kw: as_points(self, **kw),
        "as_array":       lambda self, **kw: as_array(self, **kw),
        "as_matrix":      lambda self, **kw: as_matrix(self, **kw),
        "as_data_frame":  lambda self, **kw: as_data_frame(self, **kw),
        # plot
        "plot":           lambda self, **kw: plot(self, **kw),
        "plot_rgb":       lambda self, **kw: plot_rgb(self, **kw),
    }

    # Methods to always register (Python wrappers that supersede the C++ raw method).
    _rast_force = {
        "crop", "mask", "classify", "values", "set_values",
        "aggregate", "focal", "zonal", "extract", "app", "lapp", "tapp", "sapp",
        "merge", "mosaic", "write", "plot", "plot_rgb",
        "levels", "cats", "is_factor", "as_factor",
        "is_na", "not_na", "sum", "mean", "min", "max",
    }
    for name, fn in _rast_methods.items():
        if name in _rast_force or not hasattr(SpatRaster, name):
            setattr(SpatRaster, name, fn)

    # ── SpatVector methods ────────────────────────────────────────────────────

    from .generics import (
        project_vector, shift_vect, rotate_vect, rescale_vect, trans_vect,
    )
    from .geom import (
        is_valid, make_valid,
        union_vect, intersect_vect, erase, symdif, cover_vect,
        crop_vect, mask_vect,
        buffer_vect, disagg_vect, flip_vect, spin,
        hull, delaunay, voronoi,
        simplify_geom, thin_geom, gaps,
        is_empty,
    )
    from .distance import distance_vect_self, distance_vect
    from .spatvec import (
        geomtype, is_lines, is_polygons, is_points,
        geom, crds,
        expanse, perim, nseg,
        fill_holes, vect_as_df, geom_as_wkt,
    )
    from .relate import is_related, relate, relate_self
    from .rasterize import rasterize as rasterize_fn
    from .extract import extract as extract_fn
    from .merge import merge_vect
    from .write import write_vector
    from .names import names_vect, set_names_vect

    _vect_methods = {
        "project":        lambda self, crs, **kw: project_vector(self, crs, **kw),
        "crop":           lambda self, y, **kw: crop_vect(self, y, **kw),
        "mask":           lambda self, m, **kw: mask_vect(self, m, **kw),
        "buffer":         lambda self, width, **kw: buffer_vect(self, width, **kw),
        "hull":           lambda self, **kw: hull(self, **kw),
        "voronoi":        lambda self, **kw: voronoi(self, **kw),
        "delaunay":       lambda self, **kw: delaunay(self, **kw),
        "simplify":       lambda self, tol, **kw: simplify_geom(self, tol, **kw),
        "is_valid":       lambda self: is_valid(self),
        "make_valid":     lambda self: make_valid(self),
        "is_empty":       lambda self: is_empty(self),
        "union":          lambda self, y=None, **kw: union_vect(self, y, **kw),
        "intersect":      lambda self, y, **kw: intersect_vect(self, y, **kw),
        "erase":          lambda self, y, **kw: erase(self, y, **kw),
        "disagg":         lambda self, **kw: disagg_vect(self, **kw),
        "spin":           lambda self, angle, **kw: spin(self, angle, **kw),
        "gaps":           lambda self, **kw: gaps(self, **kw),
        "expanse":        lambda self, **kw: expanse(self, **kw),
        "perim":          lambda self: perim(self),
        "geomtype":       lambda self: geomtype(self),
        "crds":           lambda self, **kw: crds(self, **kw),
        "geom":           lambda self, **kw: geom(self, **kw),
        "as_df":          lambda self: vect_as_df(self),
        "as_wkt":         lambda self: geom_as_wkt(self),
        "relate":         lambda self, y, relation, **kw: relate(self, y, relation, **kw),
        "is_related":     lambda self, y, relation, **kw: is_related(self, y, relation, **kw),
        "distance":       lambda self, **kw: distance_vect_self(self, **kw),
        "rasterize":      lambda self, r, **kw: rasterize_fn(self, r, **kw),
        "extract":        lambda self, r, **kw: extract_fn(r, self, **kw),
        "merge":          lambda self, y, **kw: merge_vect(self, y, **kw),
        "write":          lambda self, filename, **kw: write_vector(self, filename, **kw),
        "names":          lambda self: names_vect(self),
    }

    _vect_force = {
        "crop", "mask", "buffer", "simplify", "project",
        "union", "intersect", "erase", "merge", "write",
        "rasterize", "extract", "distance",
    }
    for name, fn in _vect_methods.items():
        if name in _vect_force or not hasattr(SpatVector, name):
            setattr(SpatVector, name, fn)

    # ── SpatExtent methods ────────────────────────────────────────────────────

    from .extent import ext

    _ext_methods = {
        "crop":    lambda self, y, **kw: crop(y, self, **kw),
    }

    for name, fn in _ext_methods.items():
        if not hasattr(SpatExtent, name):
            setattr(SpatExtent, name, fn)
