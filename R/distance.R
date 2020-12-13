# Author: Robert J. Hijmans
# Date : July 2019
# Version 1.0
# License GPL v3

setMethod("distance", signature(x="SpatRaster", y="missing"), 
	function(x, y, grid=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt=list())
		if (grid) {
			x@ptr <- x@ptr$gridDistance(opt)
		} else {
			x@ptr <- x@ptr$rastDistance(opt)
		}
		show_messages(x, "distance")
	}
)



setMethod("buffer", signature(x="SpatRaster"), 
	function(x, width, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt=list())
		x@ptr <- x@ptr$buffer(width, opt)
		show_messages(x, "buffer")
	}
)


setMethod("distance", signature(x="SpatRaster", y="SpatVector"), 
	function(x, y, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt=list())
		x@ptr <- x@ptr$vectDistance(y@ptr, opt)
		show_messages(x, "distance")
	}
)



setMethod("distance", signature(x="SpatVector", y="missing"), 
	function(x, y, ...) {
		nr <- nrow(x)
		ptr <- x@ptr$distance_self()
		ptr <- show_messages(ptr, "distance")
		d <- ptr$values()[[1]]
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
		ptr <- x@ptr$distance_other(y@ptr, pairwise)
		ptr <- show_messages(ptr, "distance")
		d <- ptr$values()[[1]]
		if ((nx == ny) && pairwise) {
			return(d) 
		} else {
			d <- matrix(d, nrow=nx, ncol=ny)
			return(d)
		}
	}
)


setMethod("distance", signature(x="matrix", y="matrix"), 
	function(x, y, lonlat, pairwise=FALSE, ...) {
		if (missing(lonlat)) {
			stop("you must set lonlat to TRUE or FALSE")
		}
		stopifnot(ncol(x) == 2)
		stopifnot(ncol(y) == 2)
		crs <- if (lonlat){ "+proj=longlat +datum=WGS84" } else {"+proj=utm +zone=1 +datum=WGS84"}
		x <- vect(x, crs=crs)
		y <- vect(y, crs=crs)
		distance(x, y, pairwise)
	}
)
