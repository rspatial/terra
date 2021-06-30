# Author: Robert J. Hijmans
# Date : May 2019
# Version 1.0
# License GPL v3

setMethod("ifel", signature(test="SpatRaster"), 
	function(test, yes, no, filename="", ...) {
		no_num <- FALSE
		yes_num <- FALSE
		if (!inherits(no, "SpatRaster")) {
							# logical includes default NA 
			if (!(is.numeric(no) || is.logical(no))) { 
				error("ifel", "argument 'no' must be a SpatRaster, numeric or logical")
			}
			if (length(no) > 1) warn("ifel", 'only the first element of "no" is used')
			no <- no[1]
			no_num <- TRUE
		}
		if (!inherits(yes, "SpatRaster")) {
			if (!(is.numeric(yes) || is.logical(yes))) { 
				error("ifel", "argument 'yes' must be a SpatRaster, numeric or logical")
			}
			if (length(yes) > 1) warn("ifel", 'only the first element of "yes" is used')
			yes <- yes[1]
			yes_num <- TRUE
		}
		if (no_num & yes_num) {
			stopifnot(is.numeric(yes) || is.logical(yes)) 
			return (classify(test, rbind(c(0, no[1]), c(1, yes[1]))))		
		}
		
		if (no_num) {
			no <- classify(test, rbind(c(0, no[1]), c(1, NA)))
		} else {
			no <- mask(no, test, maskvalues=TRUE)
		}
		
		if (yes_num) {
			yes <- classify(test, rbind(c(1, yes[1]), c(0, NA)))
		} else {
			yes <- mask(yes, test, maskvalues=FALSE)
		}
		
		cover(no, yes, values=NA, filename=filename, ...)
	}
)



