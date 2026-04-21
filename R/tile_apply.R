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

	# resolve cores up-front so the auto tile sizer knows about it
	ncores <- if (inherits(cores, "cluster")) length(cores)
			  else max(as.integer(cores)[1], 1L)

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

	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		cl <- cores
		doclust <- TRUE
	} else if (is.numeric(cores) && cores[1] > 1) {
		cl <- parallel::makeCluster(cores[1])
		doclust <- TRUE
		on.exit(parallel::stopCluster(cl), add=TRUE)
	} else {
		cl <- NULL
	}

	if (doclust) {
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

		# wrap x once; the workers re-open it from filenames or unpack values
		px <- wrap(x)
		# wrap any terra objects in `...` so they survive serialization
		args <- lapply(list(...), .tile_apply_wrap)

		# Build per-tile jobs: each job carries its own (outer, inner) extent
		# and output file. Workers write the cropped result to that file and retunr the path
		jobs <- mapply(function(p, f) list(outer=p$outer, inner=p$inner, file=f),
					   exts, tile_files, SIMPLIFY=FALSE)

		# inline closure so the function travels cleanly across workers
		# arg names must NOT collide with parallel::parLapply args (cl, X, fun, chunk.size)
		worker <- function(.job, .px, .ufun, .args, .wopt) {
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
			# crop to the inner (un-buffered) extent so per-tile outputs
			# do not overlap and assemble cleanly with vrt()
			if (!identical(.job$outer, .job$inner)) {
				r <- terra::crop(r, terra::ext(.job$inner), snap="near")
			}
			terra::writeRaster(r, filename=.job$file,
							   overwrite=TRUE, wopt=.wopt)
			.job$file
		}
		environment(worker) <- globalenv()

		out_files <- unlist(parallel::parLapply(cl, jobs, worker,
			.px=px, .ufun=fun, .args=args, .wopt=wopt))
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
