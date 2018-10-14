# Authors: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# Licence GPL v3



if (!isGeneric('subset')) {
	setGeneric('subset', function(x, ...)
		standardGeneric('subset')) 
}

setMethod('subset', signature(x='SpatRaster'), 
function(x, subset, filename='', overwrite=FALSE, ...) {
	if (is.character(subset)) {
		i <- stats::na.omit(match(subset, names(x)))
		if (length(i)==0) {
			stop('invalid layer names')
		} else if (length(i) < length(subset)) {
			warning('invalid layer names omitted')
		}
		subset <- i
	}
	subset <- as.integer(subset)
	if (! all(subset %in% 1:nlyr(x))) {
		stop('not a valid subset')
	}
	subset <- subset - 1;
	x@ptr <- x@ptr$subset(subset, filename, overwrite)

	if (x@ptr$warning) { warning(x@ptr$warning_message)}
	if (x@ptr$error) {	stop(x@ptr$error_message)	}			
	return(x)	
} )


setMethod("[[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})
