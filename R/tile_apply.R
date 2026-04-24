# Author: Robert J. Hijmans
# Date :  April 2026
# Version 0.1
# License GPL v3


# Coerce tile extents to a list of `list(outer=, inner=)` pairs of plain
# numeric vectors (xmin, xmax, ymin, ymax). The worker reads values on
# `outer`, then the result of `fun` is cropped to `inner` before writing,
# so the per-tile output rasters never overlap and assemble cleanly with
# vrt(). For all the explicit `tiles` shapes we trust the caller (outer ==
# inner); for the auto path we expand each inner tile by `buffer` cells,
# clamped to x's extent.
.tile_apply_extents <- function(x, tiles, cores=1, buffer=0) {

	pair <- function(o, i) list(outer = as.numeric(o), inner = as.numeric(i))
	identity_pairs <- function(m) {
		lapply(seq_len(nrow(m)), function(j) {
			v <- as.numeric(m[j, ])
			pair(v, v)
		})
	}

	if (is.null(tiles)) {
		# auto path: choose tiles, optionally with a read buffer
		inner <- getTileExtents(x, cores=cores)
		if (any(buffer > 0)) {
			outer <- getTileExtents(x, cores=cores, buffer=buffer)
			# clamp outer to the raster extent so window<- never gets a
			# read window that falls outside x.
			rext <- as.vector(ext(x))   # xmin, xmax, ymin, ymax
			out <- lapply(seq_len(nrow(inner)), function(j) {
				o <- as.numeric(outer[j, ])
				o[1] <- max(o[1], rext[1])
				o[2] <- min(o[2], rext[2])
				o[3] <- max(o[3], rext[3])
				o[4] <- min(o[4], rext[4])
				pair(o, as.numeric(inner[j, ]))
			})
			return(out)
		}
		return(identity_pairs(inner))
	}

	if (inherits(tiles, "SpatExtent")) {
		v <- as.vector(tiles)
		return(list(pair(v, v)))
	}
	if (is.list(tiles)) {
		if (length(tiles) == 0) {
			error("tile_apply", "'tiles' is empty")
		}
		out <- lapply(tiles, function(t) {
			if (inherits(t, "SpatExtent")) {
				v <- as.vector(t); pair(v, v)
			} else if (is.numeric(t) && length(t) == 4) {
				v <- as.numeric(t); pair(v, v)
			} else {
				error("tile_apply", "list elements of 'tiles' must be a SpatExtent or a length-4 numeric (xmin, xmax, ymin, ymax)")
			}
		})
		return(out)
	}
	if (is.matrix(tiles)) {
		if (ncol(tiles) != 4) {
			error("tile_apply", "matrix 'tiles' must have 4 columns: xmin, xmax, ymin, ymax")
		}
		return(identity_pairs(tiles))
	}
	if (inherits(tiles, "SpatRaster") || inherits(tiles, "SpatVector") || is.numeric(tiles)) {
		m <- getTileExtents(x, tiles)
		return(identity_pairs(m))
	}
	error("tile_apply", "'tiles' must be NULL/missing (auto from block size), a SpatExtent, a list of SpatExtents, a 4-column matrix of extents, a SpatRaster or SpatVector defining tile geometry, or a numeric (rows, cols)")
}


# wrap any terra object so that it can travel across worker boundaries.
# Recurses into plain lists so users can pass e.g. a list of SpatExtents.
.tile_apply_wrap <- function(v) {
	if (inherits(v, c("SpatRaster", "SpatVector", "SpatRasterDataset",
					  "SpatRasterCollection", "SpatVectorCollection",
					  "SpatExtent", "SpatVectorProxy", "SpatGraticule"))) {
		return(wrap(v))
	}
	if (is.list(v) && !is.object(v)) {
		return(lapply(v, .tile_apply_wrap))
	}
	v
}

.tile_apply_unwrap <- function(v) {
	if (inherits(v, "Packed")) {
		unwrap(v)
	} else if (is.list(v) && !is.object(v)) {
		lapply(v, .tile_apply_unwrap)
	} else {
		v
	}
}


# Per-tile worker for the file-backed dispatch path. Defined at file scope
# (not as a closure) so it serializes cleanly across PSOCK/forked/future
# workers without dragging tile_apply()'s local environment along. Argument
# names must NOT collide with parLapply() args (cl, X, fun, chunk.size).
# `.px` is the whole packed raster, shared across jobs; the worker windows
# it down to `.job$outer` before invoking `fun`.
.tile_apply_worker <- function(.job, .px, .ufun, .args, .wopt) {
	.unw <- function(a) {
		if (inherits(a, "Packed")) terra::unwrap(a)
		else if (is.list(a) && !is.object(a)) lapply(a, .unw)
		else a
	}
	x <- terra::unwrap(.px)
	terra::window(x) <- terra::ext(.job$outer)
	.args <- lapply(.args, .unw)
	r <- do.call(.ufun, c(list(x), .args))
	if (!inherits(r, "SpatRaster")) {
		stop("'fun' must return a SpatRaster for every tile")
	}
	# crop to the inner (un-buffered) extent so per-tile outputs do not
	# overlap and assemble cleanly with vrt()
	if (!identical(.job$outer, .job$inner)) {
		r <- terra::crop(r, terra::ext(.job$inner), snap="near")
	}
	terra::writeRaster(r, filename=.job$file, overwrite=TRUE, wopt=.wopt)
	.job$file
}


# Per-tile worker for the in-memory dispatch path. The job carries its own
# pre-cropped, packed raster slice (`.job$px`), so the worker only receives
# the pixels it needs and does not have to set a window.
.tile_apply_worker_slice <- function(.job, .ufun, .args, .wopt) {
	.unw <- function(a) {
		if (inherits(a, "Packed")) terra::unwrap(a)
		else if (is.list(a) && !is.object(a)) lapply(a, .unw)
		else a
	}
	x <- terra::unwrap(.job$px)
	.args <- lapply(.args, .unw)
	r <- do.call(.ufun, c(list(x), .args))
	if (!inherits(r, "SpatRaster")) {
		stop("'fun' must return a SpatRaster for every tile")
	}
	if (!identical(.job$outer, .job$inner)) {
		r <- terra::crop(r, terra::ext(.job$inner), snap="near")
	}
	terra::writeRaster(r, filename=.job$file, overwrite=TRUE, wopt=.wopt)
	.job$file
}


# Detect future-related `cores` values without calling the strategy: a bare
# strategy function (e.g. future::multisession) is class "future"+"function";
# a tweaked strategy (e.g. future::multisession(workers=2)) is "tweaked"+
# "future". Both should route through future.apply.
.is_future_plan <- function(cores) {
	inherits(cores, c("FutureStrategy", "tweaked")) ||
	(inherits(cores, "future") && is.function(cores))
}


# Inspect `cores` and decide which dispatch path to use. Returns a named
# list with: $kind in {"seq","cluster","future"}, $ncores (int), $cl (cluster
# or NULL), $owns_cl (TRUE iff we created the cluster and must stop it).
# The future-plan setup (installing/restoring the plan) is handled in
# tile_apply() so that on.exit is attached to the right frame; here we only
# resolve the worker count.
.tile_apply_dispatch <- function(cores) {
	if (inherits(cores, "cluster")) {
		return(list(kind="cluster", ncores=length(cores), cl=cores, owns_cl=FALSE))
	}
	if (is.character(cores) && length(cores) == 1 &&
		tolower(cores) %in% c("future", "future.apply")) {
		if (!requireNamespace("future", quietly=TRUE) ||
			!requireNamespace("future.apply", quietly=TRUE)) {
			error("tile_apply", "cores='future' requires the 'future' and 'future.apply' packages")
		}
		return(list(kind="future", ncores=max(1L, future::nbrOfWorkers()),
					cl=NULL, owns_cl=FALSE))
	}
	n <- max(as.integer(cores)[1], 1L)
	if (n > 1L) {
		return(list(kind="cluster", ncores=n, cl=NULL, owns_cl=TRUE))
	}
	list(kind="seq", ncores=1L, cl=NULL, owns_cl=FALSE)
}


tile_apply <- function(x, fun, cores=1, cpkgs=NULL, tiles=NULL, buffer=0, ...,
					   filename="", overwrite=FALSE, wopt=list(),
					   overlap_fun=NULL) {

	if (!inherits(x, "SpatRaster")) {
		error("tile_apply", "'x' must be a SpatRaster")
	}
	if (window(x)) {
		error("tile_apply", "'x' already has a window set; remove it first with 'window(x) <- NULL'")
	}
	fun <- match.fun(fun)

	# `buffer` is only honoured with auto tiles. For explicit tiles the caller
	# already controls overlap (e.g. via getTileExtents(buffer=...)).
	if (!is.null(tiles) && any(buffer > 0)) {
		warn("tile_apply", "'buffer' is only used when 'tiles' is NULL; ignoring it")
		buffer <- 0
	}

	# If the caller passed a future plan/strategy via cores=, install it
	# now and restore on exit -- so subsequent dispatch and tile sizing see
	# the new nbrOfWorkers(). After installation, the dispatch path becomes
	# the regular "future" string path.
	if (.is_future_plan(cores)) {
		if (!requireNamespace("future", quietly=TRUE) ||
			!requireNamespace("future.apply", quietly=TRUE)) {
			error("tile_apply", "a future plan was passed but 'future' and 'future.apply' are not installed")
		}
		oplan <- future::plan(cores)
		on.exit(future::plan(oplan), add=TRUE)
		cores <- "future"
	}

	# resolve dispatch + cores up-front so the auto tile sizer knows about it
	disp <- .tile_apply_dispatch(cores)
	ncores <- disp$ncores

	exts <- .tile_apply_extents(x, tiles, cores=ncores, buffer=buffer)
	ntiles <- length(exts)
	if (ntiles == 0) {
		error("tile_apply", "no tiles to process")
	}

	# Per-call temp directory for intermediate per-tile files. Always written to disk
	tdir <- tempfile(pattern="tile_apply_")
	dir.create(tdir, recursive=TRUE, showWarnings=FALSE)
	tile_files <- file.path(tdir, sprintf("tile_%05d.tif", seq_len(ntiles)))
	keep_tiles <- FALSE  # set TRUE only when the result is a VRT backed by them
	on.exit({
		if (!keep_tiles) unlink(tdir, recursive=TRUE)
	}, add=TRUE)

	# When all of x's sources are in memory, ship a per-tile slice instead
	# of the whole raster. This avoids serialising values(x) once per worker
	# (or once per chunk) and removes the need for window<- inside the
	# worker. For file-backed sources we keep the path-based contract:
	# workers parallel-read their own tile via window<-, which is strictly
	# cheaper than crop()ing on the parent and pushing pixels to workers.
	in_memory <- length(sources(x)) > 0 && all(sources(x) == "")

	build_jobs <- function() {
		if (in_memory) {
			mapply(function(p, f) {
				tile <- crop(x, ext(p$outer), snap="near")
				list(px = wrap(tile), outer = p$outer,
					 inner = p$inner, file = f)
			}, exts, tile_files, SIMPLIFY=FALSE)
		} else {
			mapply(function(p, f) list(outer=p$outer, inner=p$inner, file=f),
				   exts, tile_files, SIMPLIFY=FALSE)
		}
	}

	if (disp$kind == "cluster") {
		cl <- disp$cl
		if (is.null(cl)) {
			cl <- parallel::makeCluster(disp$ncores)
		}
		if (disp$owns_cl) {
			on.exit(parallel::stopCluster(cl), add=TRUE)
		}

		# every worker needs terra
		parallel::clusterCall(cl, function() {
			suppressPackageStartupMessages(library(terra))
			invisible(NULL)
		})
		if (!is.null(cpkgs)) {
			cpkgs <- as.character(cpkgs)
			parallel::clusterCall(cl, function(pkgs) {
				for (p in pkgs) suppressPackageStartupMessages(library(p, character.only=TRUE))
				invisible(NULL)
			}, cpkgs)
		}

		# wrap any terra objects in `...` so they survive serialization
		args <- lapply(list(...), .tile_apply_wrap)
		jobs <- build_jobs()

		if (in_memory) {
			out_files <- unlist(parallel::parLapply(cl, jobs,
				.tile_apply_worker_slice,
				.ufun=fun, .args=args, .wopt=wopt))
		} else {
			# wrap x once; the workers re-open it from filenames
			px <- wrap(x)
			out_files <- unlist(parallel::parLapply(cl, jobs,
				.tile_apply_worker,
				.px=px, .ufun=fun, .args=args, .wopt=wopt))
		}
	} else if (disp$kind == "future") {
		# Plan was already installed up-front (if the caller passed one).
		# We just dispatch to future.apply::future_lapply with the same
		# wrap/unwrap contract as the cluster path.
		args <- lapply(list(...), .tile_apply_wrap)
		jobs <- build_jobs()
		fpkgs <- unique(c("terra", as.character(cpkgs)))

		if (in_memory) {
			out_files <- unlist(future.apply::future_lapply(jobs,
				.tile_apply_worker_slice,
				.ufun=fun, .args=args, .wopt=wopt,
				future.packages=fpkgs, future.seed=TRUE))
		} else {
			px <- wrap(x)
			out_files <- unlist(future.apply::future_lapply(jobs,
				.tile_apply_worker,
				.px=px, .ufun=fun, .args=args, .wopt=wopt,
				future.packages=fpkgs, future.seed=TRUE))
		}
	} else {
		# sequential path - same disk-streaming contract as the workers
		args <- list(...)
		for (i in seq_len(ntiles)) {
			p <- exts[[i]]
			y <- deepcopy(x)
			window(y) <- ext(p$outer)
			r <- do.call(fun, c(list(y), args))
			if (!inherits(r, "SpatRaster")) {
				error("tile_apply", "'fun' must return a SpatRaster for every tile")
			}
			if (!identical(p$outer, p$inner)) {
				r <- crop(r, ext(p$inner), snap="near")
			}
			writeRaster(r, filename=tile_files[i],
						overwrite=TRUE, wopt=wopt)
			rm(r, y); gc(verbose=FALSE)
		}
		out_files <- tile_files
	}

	# Single tile shortcut: just promote the one tile file.
	if (length(out_files) == 1) {
		if (nzchar(filename)) {
			ok <- file.copy(out_files, filename, overwrite=overwrite)
			if (!isTRUE(ok)) {
				error("tile_apply", "could not copy tile output to 'filename'")
			}
			return(rast(filename))
		}
		# to survives the on.exit cleanup
		stable <- tempfile(fileext=".tif")
		file.rename(out_files, stable)
		return(rast(stable))
	}

	if (is.null(overlap_fun)) {
		# VRT for non-overlapping tiles
		vrtfile <- file.path(tdir, "all.vrt")
		v <- vrt(out_files, filename=vrtfile, overwrite=TRUE)
		if (nzchar(filename)) {
			out <- writeRaster(v, filename=filename,
							   overwrite=overwrite, wopt=wopt)
		} else {
			# tile files BACK this VRT - keep them around for the session
			keep_tiles <- TRUE
			out <- v
		}
	} else {
		rcoll <- sprc(lapply(out_files, rast))
		out <- mosaic(rcoll, fun=overlap_fun, filename=filename, overwrite=overwrite, wopt=wopt)
	}

	out
}
