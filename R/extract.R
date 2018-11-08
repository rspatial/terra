# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# Licence GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warning("cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}

setMethod("extract", signature(x="SpatRaster", y="SpatLayer"), 
function(x, y, fun="", ...) { 
    r <- x@ptr$extractLayer(y@ptr, fun)
	x <- show_messages(x, "extract")		
	r
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
function(x, i, j, ... , drop=TRUE) {
	values(x)
})


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	if (any(i > 2^.Machine$double.digits)) .big_number_warning()
	if (nargs() == 3) {
		i <- cellFromRowCol(x, i, 1:ncol(x))
	} 
	r <- x@ptr$extractCell(i)
	show_messages(x)
	matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
})

setMethod("[", c("SpatRaster", "missing", "numeric"),
function(x, i, j, ... ,drop=TRUE) {
	i <- cellFromRowCol(x, 1:nrow(x), j)
	if (any(i > 2^.Machine$double.digits)) .big_number_warning()
	r <- x@ptr$extractCell(i)
	show_messages(x)
	matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ... ,drop=TRUE) {
	i <- cellFromRowCol(x, i, j)
	if (any(i > 2^.Machine$double.digits)) .big_number_warning()
	r <- x@ptr$extractCell(i)
	show_messages(x)
	matrix(r, ncol=nlyr(x), dimnames=list(NULL, names(x)))
})

