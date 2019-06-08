# Author: Robert J. Hijmans
# Date :  January 2009
# Version 0.9
# License GPL v3


setMethod("crs", signature("SpatRaster"), 
	function(x) {
		x@ptr$crs
	}
)

setMethod("crs<-", signature("SpatRaster", "character"), 
	function(x, ..., value) {
		x@ptr$crs <- value
		return(x)
	}
)


setMethod("crs", signature("SpatRaster"), 
	function(x) {
		x@ptr$crs
	}
)


setMethod("crs<-", signature("SpatRaster", "character"), 
	function(x, ..., value) {
		x@ptr$crs <- value
		return(x)
	}
)


setMethod("crs", signature("SpatVector"), 
	function(x) {
		x@ptr$crs
	}
)

setMethod("crs<-", signature("SpatVector", "character"), 
	function(x, ..., value) {
		x@ptr$crs <- trimws(value)
		return(x)
	}
)


setMethod("isLonLat", signature("SpatRaster"), 
	function(x, ...) {
		x@ptr$isLonLat()
	}
)

setMethod("couldBeLonLat", signature("SpatRaster"), 
	function(x, ...) {
		x@ptr$couldBeLonLat()
	}
)

setMethod("isGlobalLonLat", signature("SpatRaster"), 
	function(x, ...) {
		x@ptr$isGlobalLonLat()
	}
)
