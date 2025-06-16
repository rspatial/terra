
setMethod("makeTiles", signature(x="SpatRaster"),
	function(x, y, filename="tile_.tif", extend=FALSE, na.rm=FALSE, buffer=0, overwrite=FALSE, ...) {
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
		ff
	}
)


setMethod("getTileExtents", signature(x="SpatRaster"),
	function(x, y, extend=FALSE, buffer=0) {

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

