# Author: Robert J. Hijmans
# Date : July 2019
# Version 1.0
# License GPL v3



setMethod("buffer", signature(x="SpatRaster"), 
	function(x, width, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$buffer(width, opt)
		messages(x, "buffer")
	}
)


setMethod("distance", signature(x="SpatRaster", y="missing"), 
	function(x, y, grid=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (grid) {
			x@ptr <- x@ptr$gridDistance(opt)
		} else {
			x@ptr <- x@ptr$rastDistance(opt)
		}
		messages(x, "distance")
	}
)



setMethod("distance", signature(x="SpatRaster", y="SpatVector"), 
	function(x, y, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (is.longlat(x)) {
			x@ptr <- x@ptr$vectDistanceRasterize(y@ptr, TRUE, opt)		
		} else {
			x@ptr <- x@ptr$vectDistanceDirect(y@ptr, opt)
		} 
		messages(x, "distance")
	}
)



mat2wide <- function(m, sym=TRUE, keep=NULL) {
	if (inherits(m, "dist")) { 
		# sym must be true in this case
		nr <- attr(m, "Size")
		x <- rep(1:(nr-1), (nr-1):1)
		y <- unlist(sapply(2:nr, function(i) i:nr))
		cbind(x,y, as.vector(m))
	} else {
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
}

setMethod("distance", signature(x="SpatVector", y="ANY"), 
	function(x, y, sequential=FALSE, pairs=FALSE, symmetrical=TRUE) {
		if (!missing(y)) {
			error("distance", "If 'x' is a SpatVector, 'y' should be a SpatVector or missing")
		}

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
	function(x, y, pairwise=FALSE) {
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
	function(x, y, lonlat, pairwise=FALSE) {
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
	function(x, y, lonlat, sequential=FALSE) {
		crs <- ifelse(lonlat, "+proj=longlat +datum=WGS84", 
							  "+proj=utm +zone=1 +datum=WGS84")
		x <- vect(x, crs=crs)
		distance(x, sequential)
	}
)

