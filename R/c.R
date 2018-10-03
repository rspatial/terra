# Author: Robert J. Hijmans
# Date :  October 2018
# Version 0.9
# Licence GPL v3


setMethod('c', signature(x='SpatRaster'), 
	function(x, ...) {
		for (i in list(...)) {
			if (class(i) == 'SpatRaster') {
				x@ptr <- x@ptr$addSource(i@ptr)
			}
		}
		if (x@ptr$warning) { warning(x@ptr$warning_message) }
		if (x@ptr$error) { stop(x@ptr$error_message) }
		return(x)
	}
)

