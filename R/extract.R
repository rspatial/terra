# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# Licence GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warning("cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}

setMethod("extract", signature(x="SpatRaster", y="SpatVector"), 
function(x, y, fun="", ...) { 
    r <- x@ptr$extractLayer(y@ptr, fun)
	x <- show_messages(x, "extract")		
	matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
})

setMethod("[", c("SpatRaster", "SpatVector", "missing"),
function(x, i, j, ... , drop=FALSE) {
	v <- extract(x, i)
	if (drop) {
		as.vector(v)
	} else {
		v
	}
})


setMethod("extract", signature(x="SpatRaster", y="matrix"), 
function(x, y, fun="", ...) { 
	if (ncol(y) != 2) {
		stop("extract works with a 2 column matrix of x and y coordinates")
	}
	i <- cellFromXY(x, y)
	x[i]
})


setMethod("[", c("SpatRaster", "missing", "missing"),
function(x, i, j, ... , drop=FALSE) {
	values(x, matrix=drop)
})


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	if (any(na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	if (nargs() > 2) {
		i <- cellFromRowCol(x, i, 1:ncol(x))
		# probably better to do return( readValues(x, i-1) )
	} 
	r <- x@ptr$extractCell(i-1)
	show_messages(x)
	if (drop) {
		r
	} else {
		matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
	}
})

setMethod("[", c("SpatRaster", "missing", "numeric"),
function(x, i, j, ... ,drop=FALSE) {
	i <- cellFromRowCol(x, 1:nrow(x), j)
	if (any(na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	r <- x@ptr$extractCell(i-1)
	show_messages(x)
	if (drop) {
		r
	} else {
		matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
	}
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ..., drop=FALSE) {
	i <- cellFromRowCol(x, i, j)
	if (any(na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	r <- x@ptr$extractCell(i-1)
	show_messages(x)
	if (drop) {
		r
	} else {
		matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
	}
})

