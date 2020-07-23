# Author: Robert J. Hijmans
# Date :  July 2019
# Version 1.0
# License GPL v3

setMethod("lines", signature(x="SpatExtent"), 
	function(x, col, ...)  {
		e <- as.vector(x)
		p <- rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3]), c(e[1],e[3]))
		if (missing(col)) col <- "black"
		graphics::lines(p, col=col[1], ...)		
	}
)

setMethod("plot", signature(x="SpatExtent", y="missing"), 
	function(x, y, col, ...)  {
		#e <- as.vector(x)
		#p <- rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3]), c(e[1],e[3]))
		#plot(p, type="n", axes=axes, xlab=xlab, ylab=ylab)
		#if(missing(col)) col <- "black"
		#graphics::lines(p, col=col[1], ...)		
		plot(as.polygons(x), col=col, ...)
	}
)

setMethod("points", signature(x="SpatExtent"), 
	function(x, col, ...)  {
		e <- as.vector(x)
		p <- rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3]), c(e[1],e[3]))
		if (missing(col)) col <- "black"
		col <- .getCols(4, col)
		if (is.null(col)) col <- "black"
		graphics::points(p, col=col, ...)
	}
)

