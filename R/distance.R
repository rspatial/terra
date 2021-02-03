# Author: Robert J. Hijmans
# Date : July 2019
# Version 1.0
# License GPL v3


setMethod("distance", signature(x="SpatRaster", y="missing"), 
	function(x, y, grid=FALSE, filename="", overwrite=FALSE, ...) {
		opt <- spatOptions(filename, overwrite, ...)
		if (grid) {
			x@ptr <- x@ptr$gridDistance(opt)
		} else {
			x@ptr <- x@ptr$rastDistance(opt)
		}
		messages(x, "distance")
	}
)



setMethod("buffer", signature(x="SpatRaster"), 
	function(x, width, filename="", overwrite=FALSE, ...) {
		opt <- spatOptions(filename, overwrite, ...)
		x@ptr <- x@ptr$buffer(width, opt)
		messages(x, "buffer")
	}
)


setMethod("distance", signature(x="SpatRaster", y="SpatVector"), 
	function(x, y, filename="", overwrite=FALSE, ...) {
		opt <- spatOptions(filename, overwrite, ...)
		x@ptr <- x@ptr$vectDistance(y@ptr, opt)
		messages(x, "distance")
	}
)


mat2wide <- function(m, sym=TRUE, keep=NULL) {
	bool <- is.logical(m)
	if (sym) {
		m[lower.tri(m)] <- NA
	}
	m <- cbind(from=rep(1:nrow(m), each=ncol(m)), to=rep(1:ncol(m), nrow(m)), value=as.vector(t(m)))
	m <- m[!is.na(m[,3]), , drop=FALSE]
	if (!is.null(keep)) {
		m <- m[m[,3] == keep, 1:2, drop=FALSE]
	}
	m
}


setMethod("distance", signature(x="SpatVector", y="ANY"), 
	function(x, y, sequential=FALSE, pairs=FALSE, symmetrical=TRUE, ...) {
		if (sequential) {
			return( x@ptr$distance_self(sequential))
		}
		d <- x@ptr$distance_self(sequential)
		messages(x, "distance")
		class(d) <- "dist"
		attr(d, "Size") <- nrow(x)
		attr(d, "Diag") <- FALSE
		attr(d, "Upper") <- FALSE
		attr(d, "method") <- "spatial"
		if (pairs) {
			d <- as.matrix(d)
			diag(d) <- NA
			d <- mat2wide(d, symmetrical)
		}
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
		crs <- ifelse(lonlat, "+proj=longlat +datum=WGS84", 
							  "+proj=utm +zone=1 +datum=WGS84")
		x <- vect(x, crs=crs)
		y <- vect(y, crs=crs)
		distance(x, y, pairwise)
	}
)


setMethod("distance", signature(x="matrix", y="ANY"), 
	function(x, y, lonlat, sequential=FALSE, ...) {
		crs <- ifelse(lonlat, "+proj=longlat +datum=WGS84", 
							  "+proj=utm +zone=1 +datum=WGS84")
		x <- vect(x, crs=crs)
		distance(x, sequential)
	}
)

