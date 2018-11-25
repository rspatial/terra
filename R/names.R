# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod('names', signature(x='SpatRaster'), 
	function(x) { 
		x@ptr$names
	}
)


setMethod('names<-', signature(x='SpatRaster'), 
	function(x, value)  {
		if (length(value) != nlyr(x)) {
			stop('incorrect number names')
		}
		if (! x@ptr$setNames(value)) {
			stop("cannot set these names")
		}
		if (any(names(x) != value)) {
			warning("some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)


setMethod('names', signature(x='SpatVector'), 
	function(x) { 
		x@ptr$names()
	}
)

