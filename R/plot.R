# Author: Robert J. Hijmans
# Date :  June 2008
# Version 1.0
# License GPL v3

#setMethod("plot", signature(x="SpatRaster", y="numeric"), 
#	function(x, y, maxcell=100000, xlab="", ylab="", ...)  {
#		y <- as.integer(y[1])
#		stopifnot(y>0 && y<=nlyr(x))
#		x <- x[[y]]
#		stopifnot(.hasValues(x));
#		x <- sampleRegular(x, maxcell)
#		m <- matrix(values(x), nrow=nrow(x), byrow=TRUE)	
#		require(lattice)
#		lattice::levelplot(t(m[nrow(m):1, , drop=FALSE]), xlab=xlab, ylab=ylab, ...)
#	}
#)
