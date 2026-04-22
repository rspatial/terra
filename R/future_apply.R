# Author: Robert J. Hijmans
# Date :  April 2026
# Version 0.1
# License GPL v3
#
# Thin wrappers around tile_apply() that route work through the {future}
# framework. They:
#   * tile the input raster with the same auto / explicit / buffered
#     contract as tile_apply()
#   * use the currently-active future::plan() (or one passed via `plan=`)
#   * delegate the per-tile work to terra::app() / terra::focal(), running
#     under future.apply::future_lapply() inside tile_apply()
#   * rely on wrap()/unwrap() for "filename-based shipping for large rasters
#     automatically": file-backed sources travel as paths only; large
#     in-memory rasters are written to a tempfile by wrap() before being
#     handed to the workers (see ?wrap).
#
# The user is responsible for setting up the plan, e.g.
#   future::plan(future::multisession, workers = 4)
# Sequential plans (the default) are honoured too, in which case the call
# is equivalent to a single-process tile_apply() but still goes through
# future.apply (handy for testing / progress reporters).


.future_check <- function(fname) {
	if (!requireNamespace("future", quietly=TRUE) ||
		!requireNamespace("future.apply", quietly=TRUE)) {
		error(fname, "the 'future' and 'future.apply' packages are required; install them with install.packages(c('future', 'future.apply'))")
	}
}


# Run terra::app() in tiles under a future plan.
future_app <- function(x, fun, ..., plan=NULL, cpkgs=NULL,
					   tiles=NULL, buffer=0,
					   filename="", overwrite=FALSE, wopt=list()) {
	.future_check("future_app")
	if (!inherits(x, "SpatRaster")) {
		error("future_app", "'x' must be a SpatRaster")
	}
	fun <- match.fun(fun)
	use_cores <- if (!is.null(plan)) plan else "future"
	dot_args <- list(...)

	# Per-tile work: just run app() on the windowed raster, forwarding any
	# extra args destined for the user's `fun`.
	per_tile <- function(z, ..., .ufun) terra::app(z, fun = .ufun, ...)

	args <- c(list(x,
		fun = per_tile,
		.ufun = fun,
		cores = use_cores, cpkgs = cpkgs,
		tiles = tiles, buffer = buffer,
		filename = filename, overwrite = overwrite, wopt = wopt
	), dot_args)
	do.call(tile_apply, args)
}


# Run terra::focal() in tiles under a future plan, with an automatic read
# buffer so per-tile outputs match a full-raster focal() exactly.
future_focal <- function(x, w=3, fun="sum", ..., plan=NULL, cpkgs=NULL,
						 tiles=NULL, buffer=NULL,
						 filename="", overwrite=FALSE, wopt=list()) {
	.future_check("future_focal")
	if (!inherits(x, "SpatRaster")) {
		error("future_focal", "'x' must be a SpatRaster")
	}
	use_cores <- if (!is.null(plan)) plan else "future"

	# Auto-buffer = floor(max(w)/2) cells, big enough that a per-tile focal
	# sees exactly the same neighbours as the full-raster focal would.
	if (is.null(buffer)) {
		ws <- if (is.numeric(w)) w
			  else if (is.matrix(w)) dim(w)
			  else 3L
		buffer <- max(1L, as.integer(floor(max(ws) / 2)))
	}
	# Avoid the "buffer ignored when tiles is supplied" warning when the
	# caller explicitly set tiles=...
	if (!is.null(tiles)) buffer <- 0

	per_tile <- function(z, ..., .w, .ffun) terra::focal(z, w = .w, fun = .ffun, ...)

	dot_args <- list(...)
	args <- c(list(x,
		fun = per_tile,
		.w = w, .ffun = fun,
		cores = use_cores, cpkgs = cpkgs,
		tiles = tiles, buffer = buffer,
		filename = filename, overwrite = overwrite, wopt = wopt
	), dot_args)
	do.call(tile_apply, args)
}
