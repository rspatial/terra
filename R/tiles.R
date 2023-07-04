
setMethod("makeTiles", signature(x="SpatRaster"),
	function(x, y, filename="tile_.tif", extend=FALSE, na.rm=FALSE, overwrite=FALSE, ...) {
		filename <- trimws(filename[1])
		filename <- filename[!is.na(filename)]
		if (filename == "") error("makeTiles", "filename cannot be empty")
		opt <- spatOptions(filename="", overwrite=overwrite, ...)
		if (inherits(y, "SpatRaster")) {
			ff <- x@pnt$make_tiles(y@pnt, extend[1], na.rm[1], filename, opt)
		} else if (inherits(y, "SpatVector")) {
			ff <- x@pnt$make_tiles_vect(y@pnt, extend[1], na.rm[1], filename, opt)		
		} else if (is.numeric(y)) {
			if (length(y) > 2) {
				error("makeTiles", "expected one or two numbers")
			}
			y <- rep_len(y, 2)
			y <- aggregate(rast(x), y)
			ff <- x@pnt$make_tiles(y@pnt, extend[1], na.rm[1], filename, opt)			
		} else {
			error("makeTiles", "y must be a SpatRaster or SpatVector")
		}
		messages(x, "makeTiles")
		ff
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
		} 
		r@pnt <- r@pnt$make_vrt(x, options, opt)
		messages(r, "vrt")
	}
)

