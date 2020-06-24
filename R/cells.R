# Author: Robert J. Hijmans 
# Date : June 2020
# Version 1.0
# License GPL v3

 
setMethod("cells", signature(x="SpatRaster", y="missing"), 
	function(x, y, ...) {
		# is this useful?
		which(!is.na(values(x)))
	}
)



setMethod("cells", signature("SpatRaster", "SpatVector"), 
#	function(x, y, weights=FALSE, ...) {   # in a next version
	function(x, y, ...) {
		gt <- geomtype(y)
		# will become one C++ method, but for now
		if (gt == "points") {
			g <- geom(y)
			cn <- cellFromXY(x, as.matrix(g[,c("x", "y")]))
			cbind(id=g[,1], cell=cn)
		} else {
			stop("not yet implemented")
			#x@ptr$cellFromVector(y@ptr, false) 
			#x[,2] <- x[,2] + 1
		}
	}
)

setMethod("cells", signature("SpatRaster", "SpatExtent"), 
	function(x, y, ...) {
		stop("not yet implemented")
		#x@ptr$cellFromExtent(y@ptr, false) 
		#x[,2] <- x[,2] + 1
	}
)

