
setMethod("makeTiles", signature(x="SpatRaster"), 
	function(x, y, filename="tile_.tif", extend=FALSE, na.rm=FALSE, ...) {
		filename = trimws(filename[1])
		filename <- filename[!is.na(filename)]
		if (filename == "") error("makeTiles", "filename cannot be empty")
		if (!inherits(y, "SpatRaster")) error("makeTiles", "y must be a SpatRaster")
		opt <- spatOptions(filename="", ...)
		ff <- x@ptr$make_tiles(y@ptr, extend[1], na.rm[1], filename, opt)
		messages(x, "makeTiles")
		return (ff)
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
	function(x, filename="", options=NULL, overwrite=FALSE) {
		opt <- spatOptions(filename, overwrite=overwrite)
		r <- rast()
		if (is.null(options)) {
			options=""[0]
		} else {
			if (any(substr(options, 1, 1) != "-")) {
				warn("vrt", "options that do not start with '-' are ignored")
			}
		}
		r@ptr <- r@ptr$make_vrt(x, options, opt)
		messages(r, "vrt")
	}
)

