# Author: Robert J. Hijmans
# Date :  January 2009
# Version 0.9
# License GPL v3

.proj4 <- function(x) {
	x@ptr$get_crs("prj")
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
	function(x, ...) {
		x@ptr$isLonLat()
	}
)

setMethod("isLonLat", signature("SpatVector"), 
	function(x, ...) {
		x@ptr$isLonLat()
	}
)

setMethod("couldBeLonLat", signature("SpatRaster"), 
	function(x, warn=TRUE, ...) {
		b <- x@ptr$couldBeLonLat()
		show_messages(x, "couldBeLonLat")
		b
	}
)

setMethod("couldBeLonLat", signature("SpatVector"), 
	function(x, warn=TRUE, ...) {
		b <- x@ptr$couldBeLonLat()
		show_messages(x, "couldBeLonLat")
		b
	}
)


setMethod("isGlobalLonLat", signature("SpatRaster"), 
	function(x, ...) {
		x@ptr$isGlobalLonLat()
	}
)
