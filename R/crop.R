# Author: Robert J. Hijmans
# Date : September 2018
# Version 1.0
# Licence GPL v3



setMethod('crop', signature(x='SpatRaster', y='ANY'), 
function(x, y, snap="near", filename="", overwrite=FALSE, ...) {

	if (!inherits(y, 'SpatExtent')) {
		y <- try(ext(y))
		if (class(y) == 'try-error') { stop("cannot get an extent from y") }
	}
	
	x@ptr <- x@ptr$crop(y@ptr, filename, snap, overwrite)
	if (x@ptr$warning) { warning(x@ptr$warning_message)}
	if (x@ptr$error) {	stop(x@ptr$error_message)	}			
	x	
}
)


