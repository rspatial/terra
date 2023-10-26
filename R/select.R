# Author: Robert J. Hijmans
# Date : December 2011
# Version 1.0
# License GPL v3



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
	function(x, use="rec", show=TRUE, col="cyan", draw=TRUE, ...) {
		use <- substr(tolower(use), 1, 3)
		use <- match.arg(use, c("rec", "pol"))
		scol <- ifelse(draw, "red", NA)
		if (use == "rec") {
			e <- draw(col=scol)
		#	e <- as.polygons(e)
		} else {
			e <- draw("pol", col=scol)
		}
		i <- is.related(x, e, "intersects")
		x <- x[i, ]
		if (show) {
			if (geomtype(x) == "points" || geomtype(x) == "multipoints") {
				points(x, col=col, ...)
			} else {
				lines(x, col=col, ...)
			}
		}
		x
	}
)


