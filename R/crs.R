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

setMethod("crs<-", signature("SpatRaster", "character"), 
	function(x, ..., value) {
		x@ptr$set_crs(value[1])
		show_messages(x, "crs<-")
	}
)


setMethod("crs", signature("SpatVector"), 
	function(x) {
		x@ptr$get_crs("wkt")
	}
)

setMethod("crs<-", signature("SpatVector", "character"), 
	function(x, ..., value) {
		x@ptr$set_crs(value[1])
		show_messages(x, "crs<-")
	}
)


setMethod("isLonLat", signature("SpatRaster"), 
	function(x, perhaps=FALSE, warn=TRUE, global=FALSE, ...) {
		if (perhaps) {
			a <- x@ptr$couldBeLonLat()
			if (a & global) {
				a <- x@ptr$isGlobalLonLat()
			}
			if (warn) show_messages(x, "isLonLat")
			a			
		} else if (global) {
			x@ptr$isGlobalLonLat()
		} else {
			x@ptr$isLonLat()
		}
	}
)


setMethod("isLonLat", signature("SpatVector"), 
	function(x, perhaps=FALSE, warn=TRUE, ...) {
		if (perhaps) {
			a <- x@ptr$couldBeLonLat()
			if (warn) show_messages(x, "couldBeLonLat")
			a	
		} else {
			x@ptr$isLonLat()
		}
	}
)

