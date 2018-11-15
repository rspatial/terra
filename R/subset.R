# Authors: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# Licence GPL v3



setMethod('subset', signature(x='SpatRaster'), 
function(x, subset, filename="", overwrite=FALSE, wopt=list(), ...) {
	if (is.character(subset)) {
		i <- stats::na.omit(match(subset, names(x)))
		if (length(i)==0) {
			stop('invalid layer names')
		} else if (length(i) < length(subset)) {
			warning('invalid layer names omitted')
		}
		subset <- i
	}

	subset <- as.integer(subset - 1)
	
	opt <- .runOptions(filename[1], overwrite[1], wopt)
	x@ptr <- x@ptr$subset(subset, opt)
	show_messages(x, "subset")
	return(x)	
} )


setMethod("$", "SpatRaster",  function(x, name) { subset(x, name) } )


setMethod("[[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})


setMethod("[[", c("SpatRaster", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})


