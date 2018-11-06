# Authors: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# Licence GPL v3



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
	subset <- as.integer(subset) - 1
	x@ptr <- x@ptr$subset(subset, filename, overwrite)
	.messages(x, "subset")
	return(x)	
} )


setMethod("[[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})
