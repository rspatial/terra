# Author: Robert J. Hijmans
# Date : May 2019
# Version 1.0
# Licence GPL v3

setMethod("ifel", signature(test="SpatRaster", yes="ANY", no="ANY"), 
	function(test, yes, no, filename="", ...) {
		if (!inherits(no, "SpatRaster")) {
			stopifnot(is.numeric(no))
			if (length(no) > 1) warning('only the first element of "no" is used')
			no <- reclassify(test, rbind(c(0,no[1]), c(1,NA)))	
		} else {
			no <- mask(no, test, maskvalue=TRUE)
		}
		if (!inherits(yes, "SpatRaster")) {
			stopifnot(is.numeric(yes))
			if (length(yes) > 1) warning('only the first element of "yes" is used')
			yes <- reclassify(test, rbind(c(1,yes[1]), c(0,NA)))	
		}
		cover(no, yes, filename=filename)
	}
)
