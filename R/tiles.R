
setMethod("tiles", signature(x="SpatRaster"), 
	function(x, y, filename="tile_.tif", ...) {
		stopifnot(inherits(y, "SpatRaster")) 
		y <- y[[1]]
		if (!hasValues(y)) values(y) <- 1:ncell(y)
		y <- crop(y, x, snap="out")
		v <- as.polygons(y)
		#} else if (inherits(y, "SpatVector")) {
		#	if (ncol(y) > 0) {
		#		v <- aggregate(y, colnames(y)[1])	
		#	} else {
		#		v <- y
		#		values(v) <- data.frame(tile=1:nrow(v))
		#	}
		#}
		d <- unlist(as.data.frame(v))
		filename <- filename[1]
		filename <- filename[!is.na(filename)]
		filename <- filename[filename != ""]
		if (length(filename) == 0) error("tiler", "no valid filename supplied")
		e <- paste0(".", tools::file_ext(filename))
		f <- tools::file_path_sans_ext(filename)
		ff <- paste0(f, d, e)

		for (i in 1:length(v)) {
			crop(x, v[i,], filename=ff[i], ...)
		}
		return (ff)
	}
)


setMethod("vrt", signature(x="character"), 
	function(x, filename="test.vrt", overwrite=FALSE) {
		opt <- spatOptions(filename, overwrite=overwrite)
		r <- rast()
		r@ptr <- r@ptr$make_vrt(x, opt)
		messages(r)
	}
)

