# Author: Robert J. Hijmans
# Date : June 2019
# Version 1.0
# License GPL v3

setMethod("merge", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., filename="", overwrite=FALSE, wopt=list()) { 
		opt <- .runOptions(filename, overwrite, wopt)
		dots <- list(...)
		rc <- SpatRasterCollection$new()		
		rc$add(x@ptr)
		rc$add(y@ptr)
		n <- length(dots)
		if (n > 0) {
			for (i in 1:n) {
				if (inherits(dots[i], "SpatRaster")) {
					rc$add(dots[i]@ptr)
				} else {
					name <- names(dots[i])
					cls <- class(dots[i])
					stop(paste("additional arguments should be of class 'SpatRaster'\n Found argument", name, "of class: ", cls))
				}
			}
		}
		x@ptr <- rc$merge(opt)
		show_messages(x, "merge")		
	}
)

