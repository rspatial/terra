# Author: Robert J. Hijmans
# Date :  June 2008
# Version 0.9
# License GPL v3

setMethod("readValues", signature(x="SpatRaster"), 
function(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), ...) {
	x@ptr$readValues(row-1, nrows-1, col-1, ncols-1)
}
)


setMethod("values", signature(x="SpatRaster"), 
function(x, matrix=TRUE, ...) {
	if (.hasValues(x)) {
		v <- x@ptr$getValues()
		if (matrix) {
			v <- matrix(v, ncol=nlyr(x))
			colnames(v) <- names(x)
		}
	} else {
		v <- NULL
	}
	return(v)	
}
)


setMethod("values<-", signature(x="SpatRaster", "ANY"), 
	function(x, value) {
	if (is.matrix(value)) { 
		if (nrow(value) == nrow(x)) {
			value <- as.vector(t(value))
		} else {
			value <- as.vector(value)
		} 
	} else if (is.array(value)) { 
		stopifnot(length(dim(value)) == 3)
		value <- as.vector(aperm(value, c(2,1,3)))
	}
	
	if (!(is.numeric(value) || is.integer(value) || is.logical(value))) {
		stop("value must be numeric, integer, or logical")
	}

	lv <- length(value)
	nc <- ncell(x)
	nl <- nlyr(x)
	if (lv == 1) {	
		value <- rep(value, nl * nc)
	} else {
		stopifnot((lv %% nc) == 0)
		if (lv < (nc * nl)) {
			value <- rep(value, length.out=nc*nl)
		}
	}
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
	x@ptr$filenames
}

.hasMinMax <- function(x) {
	x@ptr$hasRange
}



setMethod("sources", signature(x="SpatRaster"), 
	function(x, ...) {
		src <- x@ptr$filenames
		src[src == ""] <= "memory"
		data.frame(source=src, nlyr=x@ptr$nlyrBySource(), stringsAsFactors=FALSE)
	}
)


setMethod("minmax", signature(x="SpatRaster"), 
	function(x) {
		rmin <- x@ptr$range_min
		rmax <- x@ptr$range_max
		rbind(rmin, rmax)
	}
)


setMethod("setMinMax", signature(x="SpatRaster"), 
	function(x) {
		if (!.hasMinMax(x)) {
			x@ptr$setRange
		}
	}
)
