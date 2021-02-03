# Author: Robert J. Hijmans
# Date : May 2019
# Version 1.0
# Licence GPL v3

setMethod("ifel", signature(test="SpatRaster"), 
	function(test, yes, no, filename="", overwrite=FALSE, ...) {
		if (!inherits(no, "SpatRaster")) {
			stopifnot(is.numeric(no) || is.logical(no)) # logical includes default NA 
			if (length(no) > 1) warn("ifel", 'only the first element of "no" is used')
			no <- classify(test, rbind(c(0, no[1]), c(1, NA)))
		} else {
			no <- mask(no, test, maskvalue=TRUE)
		}
		if (!inherits(yes, "SpatRaster")) {
			stopifnot(is.numeric(yes) || is.logical(yes)) 
			if (length(yes) > 1) warn("ifel", 'only the first element of "yes" is used')
			yes <- classify(test, rbind(c(1, yes[1]), c(0, NA)))
		}
		cover(no, yes, value=NA, filename=filename, overwrite=overwrite, ...)
	}
)



