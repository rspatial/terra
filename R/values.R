# Author: Robert J. Hijmans
# Date :  June 2008
# Version 0.9
# Licence GPL v3

if (!isGeneric("setRange")) {setGeneric("setRange", function(x) standardGeneric("setRange"))}	


setMethod("values", signature(x='SpatRaster'), 
function(x, matrix=TRUE, ...) {
	if (.hasValues(x)) {
		v <- x@ptr$getValues()
		if (matrix) {
			v <- matrix(v, ncol=nlayer(x))
		}
	} else {
		v <- NULL
	}
	return(v)	
}
)


setMethod('values<-', signature(x='SpatRaster', 'numeric'), 
	function(x, value) {
	if (is.matrix(value)) { 
		if (nlayer(value) == ncol(x) && ncell(value) == nrow(x)) {
			value <- as.vector(value) 
		} else if (nrow(value) == nrow(x) && ncol(value) == ncol(x)) {
			value <- as.vector(t(value))
		} else if (ncol(value)==1 | nrow(value)==1) {
			value <- as.vector(value)
		} else {
			stop('cannot use a matrix with these dimensions')
		}
	} else if (is.array(value)) { 
		stopifnot(length(dim(value)) == 3)
		value <- as.vector(aperm(value, c(2,1,3)))
	}
	
	if (!(is.numeric(value) | is.integer(value) | is.logical(value))) {
		stop('value must be numeric, integer, or logical')
	}

	if (length(value) == 1) {	
		value <- rep(value, ncell(x))
	}

	# new pointer
	y <- rast(x)
	y@ptr$setValues(value)
	y
}
)
	

.hasValues <- function(x) {
	x@ptr$hasValues
}

.inMemory <- function(x) {
	x@ptr$inMemory
}

.filenames <- function(x) {
	x@ptr$source$filename
}

.hasRange <- function(x) {
	x@ptr$hasRange
}

setMethod('range', signature(x='SpatRaster'), 
	function(x, ..., na.rm=TRUE) {
		rmin <- x@ptr$range_min
		rmax <- x@ptr$range_max
		rbind(rmin, rmax)
	}
)


setMethod('setRange', signature(x='SpatRaster'), 
	function(x) {
		if (!.hasRange(x)) {
			x@ptr$setRange
		}
	}
)
