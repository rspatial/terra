# Author: Robert J. Hijmans
# Date : July 2019
# Version 1.0
# License GPL v3

setMethod("distance", signature(x="SpatRaster", y="missing"), 
	function(x, y, grid=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- spatOptions(filename, overwrite, wopt=list())
		if (grid) {
			x@ptr <- x@ptr$gridDistance(opt)
		} else {
			x@ptr <- x@ptr$rastDistance(opt)
		}
		messages(x, "distance")
	}
)



setMethod("buffer", signature(x="SpatRaster"), 
	function(x, width, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- spatOptions(filename, overwrite, wopt=list())
		x@ptr <- x@ptr$buffer(width, opt)
		messages(x, "buffer")
	}
)


setMethod("distance", signature(x="SpatRaster", y="SpatVector"), 
	function(x, y, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- spatOptions(filename, overwrite, wopt=list())
		x@ptr <- x@ptr$vectDistance(y@ptr, opt)
		messages(x, "distance")
	}
)



setMethod("distance", signature(x="SpatVector", y="missing"), 
	function(x, y, ...) {
		nr <- nrow(x)
		d <- x@ptr$distance_self()
		messages(x, "distance")
		class(d) <- "dist"
		attr(d, "Size") <- nr
		attr(d, "Diag") <- FALSE
		attr(d, "Upper") <- FALSE
		attr(d, "method") <- "spatial"
		d
	}
)


setMethod("distance", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, pairwise=FALSE, ...) {
		nx <- nrow(x)
		ny <- nrow(y)
		d <- x@ptr$distance_other(y@ptr, pairwise)
		messages(x, "distance")
		if ((nx == ny) && pairwise) {
			d
		} else {
			d <- matrix(d, nrow=nx, ncol=ny)
		}
		return(d)
	}
)


setMethod("distance", signature(x="matrix", y="matrix"), 
	function(x, y, lonlat, pairwise=FALSE, ...) {
		if (missing(lonlat)) {
			error("distance", "you must set lonlat to TRUE or FALSE")
		}
		stopifnot(ncol(x) == 2)
		stopifnot(ncol(y) == 2)
		crs <- if (lonlat){ "+proj=longlat +datum=WGS84" } else {"+proj=utm +zone=1 +datum=WGS84"}
		x <- vect(x, crs=crs)
		y <- vect(y, crs=crs)
		distance(x, y, pairwise)
	}
)
