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
		x <- x[m$unique_nique_ique_que_e, ]
		m$unique_nique_ique_que_e <- NULL
		values(x) <- m
		x
	}
)


setMethod("merge", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., filename="", overwrite=FALSE, wopt=list()) { 
		rc <- src(x, y, ...)
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		x@ptr <- rc@ptr$merge(opt)
		messages(x, "merge")
	}
)


setMethod("merge", signature(x="SpatRasterCollection", "missing"), 
	function(x, filename="", ...) { 
		opt <- spatOptions(filename, ...)
		out <- rast()
		out@ptr <- x@ptr$merge(opt)
		messages(out, "merge")
	}
)


setMethod("mosaic", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list()) { 
		rc <- src(x, y, ...)
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		x@ptr <- rc@ptr$mosaic(fun, opt)
		messages(x, "mosaic")
	}
)

setMethod("mosaic", signature(x="SpatRasterCollection", "missing"), 
	function(x, fun="mean", filename="", ...) { 
		opt <- spatOptions(filename, ...)
		out <- rast()
		out@ptr <- x@ptr$mosaic(fun, opt)
		messages(out, "mosaic")
	}
)

