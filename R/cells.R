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


#setMethod("cells", signature("SpatRaster", "SpatExtent"), 
#	function(x, y, ...) {
#		p <- as.polygons(y, crs=crs(x))
#		cells(x, p)[,2]
#	}
#)

setMethod("cells", signature("SpatRaster", "SpatExtent"), 
	function(x, y, ...) {
		e <- align(y, x)
		s <- res(x)/2
		e <- as.vector(y) + c(s[1], -s[1], s[2], -s[2])
		r <- rowFromY(x, e[4:3])-1
		c <- colFromX(x, e[1:2])
		cc <- c[1]:c[2]
		rr <- (r[1]:r[2]) * ncol(x)
		rep(rr, each=length(cc)) + cc
	}
)

