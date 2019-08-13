# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# License GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warning("cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}


.makeDataFrame <- function(x, v) {
	v <- data.frame(v)	
	ff <- is.factor(x)
	if (any(ff)) {
		ff <- which(ff)
		levs <- levels(x)
		for (f in ff) {
			lev <- levs[[f]]
			v[[f]] = factor(v[[f]], levels=lev$levels)
			levels(v[[f]]) = lev$labels
		}
	}
	v
}


setMethod("extract", signature(x="SpatRaster", y="SpatVector"), 
function(x, y, fun=NULL, ..., method="simple", drop=FALSE) { 
    r <- x@ptr$extractVector(y@ptr, method[1])
	x <- show_messages(x, "extract")

	#f <- function(i) if(length(i)==0) { NA } else { i }
	#r <- rapply(r, f, how="replace")

	if (!is.null(fun)) {
	  	r <- rapply(r, fun, ...)
		r <- matrix(r, nrow=nrow(y), byrow=TRUE)
		colnames(r) <- names(x)
	} else if (drop) {
		r <- unlist(r)
		r <- matrix(r, nrow=nrow(y), byrow=TRUE)
		colnames(r) <- names(x)	
	}
	r
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
function(x, y, ...) { 
	if (ncol(y) != 2) {
		stop("extract works with a 2 column matrix of x and y coordinates")
	}
	i <- cellFromXY(x, y)
	x[i]
})


setMethod("[", c("SpatRaster", "missing", "missing"),
function(x, i, j, ... , drop=FALSE) {
	values(x, mat=drop)
})


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	if (any(stats::na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	if (nargs() > 2) {
		i <- cellFromRowCol(x, i, 1:ncol(x))
		# probably better to do return( readValues(x, i-1) )
	} 
	i[i<1] <- NA
	r <- x@ptr$extractCell(i-1)
	show_messages(x)
	if (drop) {
		r
	} else {
		r <- do.call(cbind, r)
		colnames(r) = names(x)
		r
	}
})

setMethod("[", c("SpatRaster", "missing", "numeric"),
function(x, i, j, ... ,drop=FALSE) {
	i <- cellFromRowCol(x, 1:nrow(x), j)
	if (any(stats::na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	
	r <- x@ptr$extractCell(i-1)
	show_messages(x)
	if (drop) {
		r
	} else {
		r <- do.call(cbind, r)
		colnames(r) = names(x)
		r
	}
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ..., drop=FALSE) {
	i <- cellFromRowColCombine(x, i, j)
	if (any(stats::na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	r <- x@ptr$extractCell(i-1)
	show_messages(x)
	if (drop) {
		r
	} else {
		r <- do.call(cbind, r)
		colnames(r) = names(x)
		r
	}
})

