# Author: Robert J. Hijmans
# Date :  July 2019
# Version 1.0
# License GPL v3

setMethod("lines", signature(x="SpatExtent"),
	function(x, col="black", alpha=1, ...)  {
		e <- as.vector(x)
		p <- rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3]), c(e[1],e[3]))
		col <- .getCols(1, col, alpha)
		graphics::lines(p, col=col, ...)
	}
)

setMethod("polys", signature(x="SpatExtent"),
	function(x, col="black", alpha=1, ...)  {
		polys(as.polygons(x), col=col, alpha=alpha, ...)
	}
)

setMethod("plot", signature(x="SpatExtent", y="missing"),
	function(x, ...)  {
		if (!is.valid(x)) {
			error("plot", "invalid SpatExtent")
		}
		plot(as.polygons(x), ...)
	}
)

setMethod("points", signature(x="SpatExtent"),
	function(x, col="black", alpha=1, ...)  {
		e <- as.vector(x)
		p <- rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3]), c(e[1],e[3]))
		col <- .getCols(4, col, alpha)
		graphics::points(p, col=col, ...)
	}
)

