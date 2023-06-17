# Author: Robert J. Hijmans
# Date : July 2019
# Version 1.0
# License GPL v3


setMethod("buffer", signature(x="SpatRaster"),
	function(x, width, background=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$buffer(width, background, opt)
		messages(x, "buffer")
	}
)


setMethod("distance", signature(x="SpatRaster", y="missing"),
	function(x, y, target=NA, exclude=NULL, unit="m", haversine=TRUE, filename="", ...) {
		if (!is.null(list(...)$grid)) {
			error("distance", "use 'gridDistance(x)' instead of  'distance(x, grid=TRUE)'")
		}
		opt <- spatOptions(filename, ...)
		target <- as.numeric(target[1])
		keepNA <- FALSE
		if (!is.null(exclude)) {
			exclude <- as.numeric(exclude[1])
			if ((is.na(exclude) && is.na(target)) || isTRUE(exclude == target)) {
				error("distance", "'target' and 'exclude' must be different") 
			}
			if (is.na(exclude)) {
				keepNA <- TRUE
			}
		} else {
			exclude <- NA
		}
		x@ptr <- x@ptr$rastDistance(target, exclude, keepNA, tolower(unit), TRUE, haversine, opt)
		messages(x, "distance")
	}
)


setMethod("costDist", signature(x="SpatRaster"),
	function(x, target=0, scale=1, maxiter=50, filename="", ...) {
		opt <- spatOptions(filename, ...)
		maxiter <- max(maxiter[1], 2)
		x@ptr <- x@ptr$costDistance(target[1], scale[1], maxiter, FALSE, opt)
		messages(x, "costDist")
	}
)


setMethod("gridDist", signature(x="SpatRaster"),
	function(x, target=0, scale=1, maxiter=50, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (is.na(target)) {
			x@ptr <- x@ptr$gridDistance(scale[1]	, opt)
		} else {
			maxiter <- max(maxiter[1], 2)
			x@ptr <- x@ptr$costDistance(target[1], scale[1], maxiter, TRUE, opt)
		}
		messages(x, "gridDist")
	}
)


setMethod("distance", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, unit="m", rasterize=FALSE, haversine=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		unit <- as.character(unit[1])
		if (rasterize) {
			x@ptr <- x@ptr$vectDistanceRasterize(y@ptr, NA, NA, unit, haversine, opt)
		} else {
			x@ptr <- x@ptr$vectDistanceDirect(y@ptr, unit, haversine, opt)
		}
		messages(x, "distance")
	}
)


setMethod("distance", signature(x="SpatRaster", y="sf"),
	function(x, y, unit="m", rasterize=FALSE, haversine=TRUE, filename="", ...) {
		distance(x, vect(y), unit=unit, rasterize=rasterize, haversine=haversine, filename=filename, ...) 
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
	function(x, y, sequential=FALSE, pairs=FALSE, symmetrical=TRUE, unit="m") {
		if (!missing(y)) {
			error("distance", "If 'x' is a SpatVector, 'y' should be a SpatVector or missing")
		}

		if (sequential) {
			return( x@ptr$distance_self(sequential, unit))
		}
		unit <- as.character(unit[1])
		d <- x@ptr$distance_self(sequential, unit)
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
	function(x, y, pairwise=FALSE, unit="m") {
		unit <- as.character(unit[1])
		d <- x@ptr$distance_other(y@ptr, pairwise, unit)
		messages(x, "distance")
		if (!pairwise) {
			d <- matrix(d, nrow=nrow(x), ncol=nrow(y), byrow=TRUE)
		}
		d
	}
)


setMethod("distance", signature(x="matrix", y="matrix"),
	function(x, y, lonlat, pairwise=FALSE) {
		if (missing(lonlat)) {
			error("distance", "you must set lonlat to TRUE or FALSE")
		}
		stopifnot(ncol(x) == 2)
		stopifnot(ncol(y) == 2)
		v <- vect()
		d <- v@ptr$point_distance(x[,1], x[,2], y[,1], y[,2], pairwise[1], 1, lonlat)
		messages(v)
		if (pairwise) {
			d
		} else {
			matrix(d, nrow=nrow(x), ncol=nrow(y), byrow=TRUE)
		}
	}
)


setMethod("distance", signature(x="matrix", y="missing"),
	function(x, y, lonlat=NULL, sequential=FALSE, pairs=FALSE, symmetrical=TRUE) {
		if (is.null(lonlat)) {
			error("distance", "lonlat should be TRUE or FALSE")
		}
		crs <- ifelse(isTRUE(lonlat), "+proj=longlat +datum=WGS84",
							          "+proj=utm +zone=1 +datum=WGS84")
		x <- vect(x, crs=crs)
		distance(x, sequential=sequential, pairs=pairs, symmetrical=symmetrical)
	}
)


setMethod("direction", signature(x="SpatRaster"),
	function(x, from=FALSE, degrees=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rastDirection(from[1], degrees[1], NA, NA, opt)
		messages(x, "direction")
	}
)

