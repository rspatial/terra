# Author: Robert J. Hijmans
# Date :  June 2008
# Version 1.0
# Licence GPL v3


if (!isGeneric("plot")) {
	setGeneric("plot", function(x,y,...)
		standardGeneric("plot"))
}	

setMethod("plot", signature(x='SpatRaster', y='missing'), 
	function(x, y, maxpixels=500000, xlab="", ylab="", ...)  {
		stopifnot(.hasValues(x));
		require(lattice)
		m <- values(x)[, 1]		
		m <- matrix(m, nrow=nrow(x), byrow=TRUE)
		lattice::levelplot(t(m[nrow(m):1, , drop=FALSE]), xlab=xlab, ylab=ylab, ...)
	}
)


setMethod("plot", signature(x='SpatRaster', y='numeric'), 
	function(x, y, maxpixels=500000, xlab="", ylab="", ...)  {
		stopifnot(.hasValues(x));
		require(lattice)
		m <- values(x)[, y]
		m <- matrix(m, nrow=nrow(x), byrow=TRUE)
		lattice::levelplot(t(m[nrow(m):1, , drop=FALSE]), xlab=xlab, ylab=ylab, ...)
	}
)


