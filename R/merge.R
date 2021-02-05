# Author: Robert J. Hijmans
# Date : June 2019
# Version 1.0
# License GPL v3


setMethod("mosaic", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list()) { 
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		dots <- list(...)
		rc <- SpatRasterCollection$new()
		rc$add(x@ptr)
		rc$add(y@ptr)
		n <- length(dots)
		if (n > 0) {
			for (i in 1:n) {
				if (inherits(dots[[i]], "SpatRaster")) {
					rc$add(dots[[i]]@ptr)
				} else {
					name <- names(dots[[i]])
					cls <- class(dots[[i]])
					error("mosaic", "additional arguments should be 'SpatRaster'\n Found argument", name, "of class: ", cls)
				}
			}
		}
		x@ptr <- rc$mosaic(fun, opt)
		messages(x, "mosaic")
	}
)

setMethod("merge", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., filename="", overwrite=FALSE, wopt=list()) { 
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		dots <- list(...)
		rc <- SpatRasterCollection$new()
		rc$add(x@ptr)
		rc$add(y@ptr)
		n <- length(dots)
		if (n > 0) {
			for (i in 1:n) {
				if (inherits(dots[[i]], "SpatRaster")) {
					rc$add(dots[[i]]@ptr)
				} else {
					name <- names(dots[[i]])
					cls <- class(dots[[i]])
					error("merge", "additional arguments should be 'SpatRaster'\n Found argument", name, "of class: ", cls)
				}
			}
		}
		x@ptr <- rc$merge(opt)
		messages(x, "merge")
	}
)


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



