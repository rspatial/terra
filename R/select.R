# Author: Robert J. Hijmans
# Date : December 2011
# Version 1.0
# Licence GPL v3



setMethod("sel", signature(x="SpatRaster"), 
	function(x, ...) {
		e <- draw(...)
		int <- intersect(e, ext(x))
		if (is.null(int)) {
			x <- NULL
		} else {
			x <- crop(x, e)
		}
		x
	}
)


setMethod("sel", signature(x="SpatVector"), 
	function(x, use="rec", draw=TRUE, col="cyan", ...) {
		use <- substr(tolower(use), 1, 3)
		use <- match.arg(use, c("rec", "pol")) 
		if (use == "rec") {
			e <- draw()
		#	e <- as.polygons(e)
		} else {
			e <- draw("pol")
		}
		i <- relate(x, e, "intersects")
		x <- x[as.vector(i), ]
		if (draw) {
			if (geomtype(x) == "points" || geomtype(x) == "multipoints") {
				points(x, col=col, ...)
			} else {
				lines(x, col=col, ...)
			}
		}
		x
	}
)


