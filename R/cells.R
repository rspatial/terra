# Author: Robert J. Hijmans
# Date : June 2020
# Version 1.0
# License GPL v3


setMethod("cells", signature(x="SpatRaster", y="missing"),
	function(x, y) {
		if (hasValues(x)) {
			opt <- spatOptions()
			x@ptr$cells_notna_novalues(opt) + 1
		} else {
			# is this useful?
			1:ncell(x)
		}
	}
)

setMethod("cells", signature(x="SpatRaster", y="numeric"),
	function(x, y, pairs=FALSE) {
		opt <- spatOptions()
		v <- x@ptr$is_in_cells(y, pairs, opt)
		x <- messages(x, "cells")
		if (pairs) {
			v <- lapply(v, function(i) {
				m <- matrix(i, ncol=2)
				m[,1] <- m[,1] + 1
				colnames(m) <- c("cell", "value")
				m
			})
		} else {
			v <- lapply(v, function(i) i+1)
		}
		names(v) <- names(x)
		v
	}
)


setMethod("cells", signature("SpatRaster", "SpatVector"),
	function(x, y, method="simple", weights=FALSE, exact=FALSE, touches=is.lines(y), small=TRUE) {
		method = match.arg(tolower(method), c("simple", "bilinear"))
		opt <- spatOptions()
		d <- x@ptr$vectCells(y@ptr, touches[1], small[1], method[1], weights[1], exact[1], opt)
		if (geomtype(y) == "points") {
			d <- matrix(d, nrow=nrow(y), byrow=TRUE)
			d <- cbind(1:nrow(y), d)
			if (method == "bilinear") {
				colnames(d) <- c("ID", "c1", "c2", "c3", "c4", "w1", "w2", "w3", "w4")
				d[,2:5] <- d[,2:5] + 1
			} else {
				colnames(d) <- c("ID", "cell")
				d[,2] <- d[,2] + 1
			}
			return (d)
		}
		cn <- c("ID", "cell")
		if (weights[1] || exact[1]) {
			d <- matrix(d, ncol=3)
			cn <- c(cn, "weights")
		} else {
			d <- matrix(d, ncol=2)
		}
		d[,1:2] <- d[,1:2] + 1
		colnames(d) <- cn
		d
	}
)


#setMethod("cells", signature("SpatRaster", "SpatExtent"),
#	function(x, y, ...) {
#		p <- as.polygons(y, crs=crs(x))
#		cells(x, p)[,2]
#	}
#)

#setMethod("cells", signature("SpatRaster", "SpatExtent"),
#	function(x, y, ...) {
#		e <- align(y, x)
#		s <- res(x)/2
#		e <- as.vector(y) + c(s[1], -s[1], s[2], -s[2])
#		r <- rowFromY(x, e[4:3])-1
#		c <- colFromX(x, e[1:2])
#		cc <- c[1]:c[2]
#		rr <- (r[1]:r[2]) * ncol(x)
#		rep(rr, each=length(cc)) + cc
#	}
#)


setMethod("cells", signature("SpatRaster", "SpatExtent"),
	function(x, y) {
		opt <- spatOptions()
		x@ptr$extCells(y@ptr) + 1
	}
)
