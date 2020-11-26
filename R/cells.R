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

setMethod("cells", signature(x="SpatRaster", y="numeric"), 
	function(x, y, ...) {
		opt <- .runOptions("", TRUE, list())
		v <- x@ptr$is_in_cells(y, opt)
		x <- show_messages(x, "cells")
		v <- lapply(v, function(i) i+1)
		names(v) <- names(x)
		v
	}
)


setMethod("cells", signature("SpatRaster", "SpatVector"), 
#	function(x, y, weights=FALSE, touches=is.lines(y), method="simple", ...) {   # in a next version
	function(x, y, touches=is.lines(y), method="simple", ...) {
		d <- x@ptr$getCells(y@ptr, touches[1], method) 
		d <- matrix(d, ncol=2)
		colnames(d) <- c("id", "cell")
		d[,2] <- d[,2] + 1
		d
	}
)


setMethod("cells", signature("SpatRaster", "SpatExtent"), 
	function(x, y, ...) {
		p <- as.polygons(y)
		cells(x, p)[,2]
	}
)

