# Author: Robert J. Hijmans
# Date : June 2019
# Version 1.0
# License GPL v3

setMethod("merge", signature(x="SpatVector", y="data.frame"),
	function(x, y, ...) {
		v <- values(x)
		v$unique_nique_ique_que_e <- 1:nrow(v)
		m <- merge(v, y, ...)
		m <- m[order(m$unique_nique_ique_que_e), ]
		x <- x[stats::na.omit(m$unique_nique_ique_que_e), ]
		m$unique_nique_ique_que_e <- NULL
		if (nrow(m) > nrow(x)) {
			error("merge", "using 'all.y=TRUE' is not allowed. Should it be?")
		}
		values(x) <- m
		x
	}
)

setMethod("merge", signature(x="SpatVector", y="SpatVector"),
	function(x, y, ...) {
		merge(x, data.frame(y), ...)
	}
)

setMethod("merge", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, ..., first=TRUE, na.rm=TRUE, filename="", overwrite=FALSE, wopt=list()) {
		rc <- sprc(x, y, ...)
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		x@cpp <- rc@cpp$merge(first[1], na.rm, opt)
		messages(x, "merge")
	}
)


setMethod("merge", signature(x="SpatRasterCollection", "missing"),
	function(x, first=TRUE, na.rm=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		out <- rast()
		out@cpp <- x@cpp$merge(first[1], na.rm, opt)
		messages(out, "merge")
	}
)


setMethod("mosaic", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list()) {
		fun <- .makeTextFun(fun)
		if (!inherits(fun, "character")) {
			error("mosaic", "function 'fun' is not valid")
		}
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		rc <- sprc(x, y, ...)
		x@cpp <- rc@cpp$mosaic(fun, opt)
		messages(x, "mosaic")
	}
)

setMethod("mosaic", signature(x="SpatRasterCollection", "missing"),
	function(x, fun="mean", filename="", ...) {
		opt <- spatOptions(filename, ...)
		out <- rast()
		fun <- .makeTextFun(fun)
		if (!inherits(fun, "character")) {
			error("mosaic", "function 'fun' is not valid")
		}
		out@cpp <- x@cpp$mosaic(fun, opt)
		messages(out, "mosaic")
	}
)

