

setMethod("makeTiles", signature(x="SpatRaster"),
	function(x, y, filename="tile_.tif", extend=FALSE, na.rm=FALSE, buffer=0, value="files", overwrite=FALSE, ...) {
	
		value <- match.arg(value, c("files", "raster", "collection"))
		filename <- trimws(filename[1])
		filename <- filename[!is.na(filename)]
		if (filename == "") error("makeTiles", "filename cannot be empty")
		opt <- spatOptions(filename="", overwrite=overwrite, ...)
		if (inherits(y, "SpatRaster")) {
			ff <- x@pntr$make_tiles(y@pntr, extend[1], buffer, na.rm[1], filename, opt)
		} else if (inherits(y, "SpatVector")) {
			ff <- x@pntr$make_tiles_vect(y@pntr, extend[1], buffer, na.rm[1], filename, opt)
		} else if (is.numeric(y)) {
			if (length(y) > 2) {
				error("makeTiles", "expected one or two numbers")
			}
			y <- rep_len(y, 2)
			y <- aggregate(rast(x), y)
			ff <- x@pntr$make_tiles(y@pntr, extend[1], buffer, na.rm[1], filename, opt)
		} else {
			error("makeTiles", "y must be numeric or a SpatRaster or SpatVector")
		}
		messages(x, "makeTiles")
		if (value == "files") {
			ff
		} else if (value == "raster") {
			filename <- gsub("_\\.", "s.", filename)
			vrt(ff, filename=filename)
		} else {
			sprc(ff)
		}
	}
)


# Choose a sensible (rows, cols) tile size for x.
#  - If x reports a GDAL block size for its source(s) (typical for tiled
#    GeoTIFFs, COGs, NetCDF, ...), align tiles to whole blocks; otherwise
#    fall back to a 256-cell square (a common GDAL block default).
#  - Cap the per-tile cell count by a memory budget, so that `cores`
#    workers running concurrently fit in `memfrac` of free RAM.
#  - Aim for a few tasks per worker (`tasks_per_worker`) to load-balance.
#  - Result is clamped to the raster's own dimensions.
.auto_tile_size <- function(x, cores=1, memfrac=NULL,
							tasks_per_worker=4, ncopies=4) {

	nr <- nrow(x); nc <- ncol(x); nl <- nlyr(x)
	cores <- max(as.integer(cores)[1], 1L)
	tasks_per_worker <- max(as.integer(tasks_per_worker)[1], 1L)

	# block size from GDAL (per source). Use the smallest reported per dim,
	# then fall back to 256 when nothing useful is reported (e.g. for
	# in-memory rasters or sources that report 0).
	fb <- tryCatch(fileBlocksize(x), error=function(e) NULL)
	if (is.null(fb) || NROW(fb) == 0) {
		br <- bc <- 0L
	} else {
		br <- as.integer(min(fb[, "rows"]))
		bc <- as.integer(min(fb[, "cols"]))
	}
	if (is.na(br) || br <= 0) br <- min(256L, nr)
	if (is.na(bc) || bc <= 0) bc <- min(256L, nc)
	br <- min(br, nr); bc <- min(bc, nc)
	one_block <- as.numeric(br) * bc

	# memory budget (bytes per cell assumes double precision and `ncopies`
	# in-flight copies, the same accounting blocks() and writeStart() use).
	bytes_per_cell <- 8 * max(nl, 1L) * max(ncopies, 1L)

	if (is.null(memfrac)) {
		memfrac <- tryCatch(spatOptions()$memfrac, error=function(e) 0.5)
		if (is.null(memfrac) || !is.finite(memfrac) || memfrac <= 0) {
			memfrac <- 0.5
		}
	}

	free_kb <- tryCatch(free_RAM(), error=function(e) NA_real_)
	if (is.na(free_kb) || !is.finite(free_kb) || free_kb <= 0) {
		# no memory info: only the parallelism constraint applies
		mem_cells <- Inf
	} else {
		free_b <- free_kb * 1024              # free_RAM returns kB
		budget <- (free_b * memfrac) / cores  # per-worker peak budget
		mem_cells <- budget / bytes_per_cell
	}

	# at least one tile per worker, and ideally a few for load balancing.
	# Avoid integer overflow on large rasters by going through `numeric`.
	total_cells <- as.numeric(nr) * nc
	target_min_tiles <- cores * tasks_per_worker
	par_cells <- total_cells / target_min_tiles

	cells_per_tile <- min(mem_cells, par_cells)
	# never go below a single block (avoids dicing tiled COGs into sub-blocks)
	cells_per_tile <- max(cells_per_tile, one_block)
	# never go above the whole raster
	cells_per_tile <- min(cells_per_tile, total_cells)

	# pick (rows, cols) as integer multiples of (br, bc), as square as possible
	side <- sqrt(cells_per_tile)
	nbr <- max(round(side / br), 1)
	nbc <- max(round(side / bc), 1)
	tile_rows <- min(as.integer(nbr * br), nr)
	tile_cols <- min(as.integer(nbc * bc), nc)

	c(tile_rows, tile_cols)
}


setMethod("getTileExtents", signature(x="SpatRaster"),
	function(x, y, extend=FALSE, buffer=0, cores=1) {

		if (missing(y) || is.null(y)) {
			y <- .auto_tile_size(x, cores=cores)
		}

		opt <- spatOptions(filename="")
		if (inherits(y, "SpatRaster")) {
			e <- x@pntr$get_tiles_ext(y@pntr, extend[1], buffer)
		} else if (inherits(y, "SpatVector")) {
			e <- x@pntr$get_tiles_ext_vect(y@pntr, extend[1], buffer)
		} else if (is.numeric(y)) {
			if (length(y) > 2) {
				error("getTileExtents", "expected one or two numbers")
			}
			y <- rep_len(y, 2)
			y <- aggregate(rast(x), y)
			e <- x@pntr$get_tiles_ext(y@pntr, extend[1], buffer)
		} else {
			error("getTileExtents", "y must be numeric or a SpatRaster or SpatVector")
		}
		messages(x, "getTileExtents")
		e <- matrix(e, ncol=4, byrow=FALSE)
		colnames(e) <- c("xmin", "xmax", "ymin", "ymax")
		e
	}
)



#		if (!hasValues(x)) error("makeTiles", "x has no values")
#		y <- rast(y)[[1]]
#		if (expand) y <- expand(y, ext(x), snap="out")
#		y <- crop(rast(y)[[1]], x, snap="out")
#		d <- 1:ncell(y)
#		if (length(filename) == 0) error("tiler", "no valid filename supplied")
#		e <- paste0(".", tools::file_ext(filename))
#		f <- tools::file_path_sans_ext(filename)
#		ff <- paste0(f, d, e)
#		for (i in d) {
#			crop(x, y[i,drop=FALSE], filename=ff[i], ...)
#		}
#		ff[file.exists(ff)]
#	}
#)


setMethod("vrt", signature(x="character"),
	function(x, filename="", options=NULL, overwrite=FALSE, set_names=FALSE, return_filename=FALSE) {
		opt <- spatOptions(filename, overwrite=overwrite)
		r <- rast()
		if (is.null(options)) {
			options=""[0]
		}
		f <- r@pntr$make_vrt(x, options, opt)
		messages(r, "vrt")
		messages(opt, "vrt")
		if (set_names) {
			v <- readLines(f)
			if ("-separate" %in% options) {
			    nms <- unlist(lapply(sprc(x), names))
			} else {
			    nms <- names(rast(x[1]))
			}
			i <- grep("band=", v)
			if (length(i) == length(nms)) {
				nms <- paste0("<Description>", nms, "</Description>")
				v[i] <- paste(v[i], nms)
				writeLines(v, f)
			}
		}
		if (return_filename) {
			f
		} else {
			rast(f)
		}
	}
)


setMethod("vrt", signature(x="SpatRasterCollection"),
	function(x, filename="", options=NULL, overwrite=FALSE, return_filename=FALSE) {
		opt <- spatOptions(filename, overwrite=overwrite)
		if (is.null(options)) {
			options=""[0]
		}
		f <- x@pntr$make_vrt(options, FALSE, opt)
		messages(x, "vrt")
		if (return_filename) {
			f
		} else {
			rast(f)
		}
	}
)


vrt_tiles <- function(x) {
	if (inherits(x, "SpatRaster")) {
		x <- sources(x)
	}
	if (!inherits(x, "character")) {
		error("vrt_sources", "x must be a filename (character) or SpatRaster)")
	}
	x <- grep(".vrt$", x, ignore.case =TRUE, value=TRUE)
	if (length(x) == 0) {
		error("vrt_sources", 'no filenames with extension ".vrt"')
	}
	tiles <- lapply(x, function(f) {
			v <- readLines(f)
			v <- v[grep("SourceFilename", v)]
			s <- strsplit(v, "\"")
			rel <- sapply(s, function(x) x[2])
			ff <- strsplit(sapply(s, function(x) x[3]), "<")
			ff <- gsub(">", "", sapply(ff, function(x) x[1]))
			ff[rel=="1"] <- file.path(dirname(f), ff[rel=="1"])
			ff
		})
	unlist(tiles)
}

