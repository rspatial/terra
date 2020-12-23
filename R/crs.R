# Author: Robert J. Hijmans
# Date :  January 2009
# Version 0.9
# License GPL v3

.proj4 <- function(x) {
	x@ptr$get_crs("proj4")
}


setMethod("crs", signature("SpatRaster"), 
	function(x) {
		x@ptr$get_crs("wkt")
	}
)


.txtCRS <- function(x, warn=TRUE) {
	if (is.na(x)) {
		x <- ""
	} else if (inherits(x, "CRS")) {
		if (warn) warn("crs", "expected a character string, not a CRS object")
		y <- attr(x, "comment")
		if (is.null(y)) {
			y <- x@projargs				
			if (is.na(y)) y <- ""
		}
		x <- y
	} else if (is.character(x)) {
		x <- x[1]
	} else {
		error("crs", "I do not know what to do with this argument (expected a character string)")
	}
	x
}

setMethod("crs<-", signature("SpatRaster", "ANY"), 
	function(x, ..., value) {
		value <- .txtCRS(value)
		x@ptr$set_crs(value)
		messages(x, "crs<-")
	}
)

#setMethod("crs<-", signature("SpatRaster", "character"), 
#	function(x, ..., value) {
#		x@ptr$set_crs(value[1])
#		messages(x, "crs<-")
#	}
#)


setMethod("crs", signature("SpatVector"), 
	function(x) {
		x@ptr$get_crs("wkt")
	}
)

setMethod("crs<-", signature("SpatVector", "ANY"), 
	function(x, ..., value) {
		value <- .txtCRS(value)
		x@ptr$set_crs(value[1])
		messages(x, "crs<-")
	}
)



setMethod("is.lonlat", signature("SpatRaster"), 
	function(x, perhaps=FALSE, warn=TRUE, global=FALSE, ...) {
		if (perhaps) {
			ok <- x@ptr$isGeographic()
			if (ok) {
				if (global) {
					return(x@ptr$isGlobalLonLat())
				} else {
					return(ok)
				}
			}
			ok <- x@ptr$couldBeLonLat()
			if (ok) {
				if (global) {
					ok <- x@ptr$isGlobalLonLat()
				}
			}
			if (ok && warn) {
				warn("is.lonlat", "assuming lon/lat crs")
			}
			return(ok)	
		} else {
			ok <- x@ptr$isGeographic()
			if (ok && global) {
				ok <- x@ptr$isGlobalLonLat()
			}
			return(ok)
		}
	}
)


setMethod("is.lonlat", signature("SpatVector"), 
	function(x, perhaps=FALSE, warn=TRUE, ...) {
		ok <- x@ptr$isGeographic()
		if (ok) return(ok)
		if (perhaps) {
			ok <- x@ptr$couldBeLonLat()
			if (ok && warn) {
				if (crs(x) == "") warn("is.lonlat", "assuming lon/lat crs")
			}
			ok
		}
	}
)

