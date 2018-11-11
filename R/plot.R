# Author: Robert J. Hijmans
# Date :  June 2008
# Version 1.0
# Licence GPL v3



setMethod("plot", signature(x='SpatRaster', y='numeric'), 
	function(x, y, maxpixels=100000, xlab="", ylab="", ...)  {
		y <- as.integer(y[1])
		stopifnot(y>0 && y<=nlyr(x))
		x <- x[[y]]
		stopifnot(.hasValues(x));
		x <- sampleRegular(x, maxpixels)
		m <- matrix(values(x), nrow=nrow(x), byrow=TRUE)	
		require(lattice)
		lattice::levelplot(t(m[nrow(m):1, , drop=FALSE]), xlab=xlab, ylab=ylab, ...)
	}
)

setMethod("plot", signature(x='SpatRaster', y='missing'), 
	function(x, y, maxpixels=100000, xlab="", ylab="", ...)  {
		plot(x, 1, maxpixels=maxpixels, xlab=xlab, ylab=ylab, ...)
	}
)

setMethod("plot", signature(x='SpatLayer', y='missing'), 
	function(x, y, xlab="", ylab="", ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		if (gtype == "points") {
			plot(g[,3], g[,4], xlab=xlab, ylab=ylab, ...)
		} else {
			g <- split(g, g[,1])
			g <- lapply(g, function(x) split(x, x[,2]))
			e <- matrix(as.vector(ext(x)), 2)
			plot(e, type='n', xlab=xlab, ylab=ylab, ...)
			p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4])))
		}
	}
)


setMethod("lines", signature(x='SpatLayer'), 
	function(x, ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		if (gtype == "points") {
			points(g[,3:4], ...)
		} else {
			g <- split(g, g[,1])
			g <- lapply(g, function(x) split(x, x[,2]))
			e <- matrix(as.vector(ext(x)), 2)
			p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4])))
		}
	}
)

setMethod("points", signature(x='SpatLayer'), 
	function(x, ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		points(g[,3:4], ...)
	}
)

