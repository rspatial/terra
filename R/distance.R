# Author: Robert J. Hijmans
# Date : July 2019
# Version 1.0
# License GPL v3

setMethod("distance", signature(x="SpatRaster", y="missing"), 
	function(x, y, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt=list())
		x@ptr <- x@ptr$rastDistance(opt)
		show_messages(x, "distance")
	}
)


setMethod("distance", signature(x="SpatRaster", y="SpatVector"), 
	function(x, y, wopt=list(), filename="", overwrite=FALSE, ...) {
		opt <- .runOptions(filename, overwrite, wopt=list())
		x@ptr <- x@ptr$pointDistance(y@ptr, opt)
		show_messages(x, "distance")
	}
)

