# Author: Robert J. Hijmans
# Date : December 2011
# Version 1.0
# Licence GPL v3



setMethod("select", signature(x="SpatRaster"), 
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
	
	
#setMethod("select", signature(x="SpatVector"), 
#	function(x, use="rec", draw=TRUE, col="cyan", size=2, ...) {
#		use <- substr(tolower(use), 1, 3)
#		stopifnot(use %in% c("rec", "pol"))
#		if (use == "rec") {
#			e <- as(drawExtent(), "SpatialPolygons")
#		} else {
#			e <- draw("pol")
#		}
#		intersect(x, e)
#	}
#)


