# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("names", signature(x="SpatRaster"), 
	function(x) { 
		x@ptr$names
	}
)


setMethod("names<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- as.character(value)
		if (length(value) != nlyr(x)) {
			stop("incorrect number of names")
		}
		if (! x@ptr$setNames(value, FALSE)) {
			stop("cannot set these names")
		}
		
		if (any(names(x) != value)) {
			# should only be possible with $setNames(value, TRUE)
			warning("some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)

setMethod("names", signature(x="SpatDataSet"), 
	function(x) { 
		x@ptr$names
	}
)


setMethod("names<-", signature(x="SpatDataSet"), 
	function(x, value) {
		value <- as.character(value)
		x@ptr$names <- value
		x
	}
)

setMethod("names", signature(x="SpatVector"), 
	function(x) { 
		x@ptr$names
	}
)


setMethod("names<-", signature(x="SpatVector"), 
	function(x, value)  {
		if (length(value) != ncol(x)) {
			stop("incorrect number names")
		}
		x@ptr$names <- value
		if (any(names(x) != value)) {
			warning("some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)
