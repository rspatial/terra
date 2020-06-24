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
#	function(x, y, weights=FALSE, touches=false, ...) {   # in a next version
	function(x, y, touches=is.lines(y), method="simple", ...) {
		d <- x@ptr$getCells(y@ptr, touches[1], method) 
		d <- matrix(d + 1, ncol=2)
		colnames(d) <- c("id", "cell")
		d
	}
)

setMethod("cells", signature("SpatRaster", "SpatExtent"), 
	function(x, y, ...) {
		p <- as.polygons(y)
		cells(x, p)[,2]
	}
)

