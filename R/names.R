# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# Licence GPL v3


setMethod('names', signature(x='SpatRaster'), 
	function(x) { 
		x@ptr$names
	}
)


setMethod('names<-', signature(x='SpatRaster'), 
	function(x, value)  {
		nl <- nlyr(x)
		if (length(value) != nl) {
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

