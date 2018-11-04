# Author: Robert J. Hijmans
# Date :  January 2009
# Version 0.9
# Licence GPL v3


setMethod("crs", signature('SpatRaster'), 
	function(x) {
		x@ptr$crs
	}
)

setMethod("crs<-", signature('SpatRaster', 'character'), 
	function(x, ..., value) {
		x@ptr$crs <- value
		return(x)
	}
)


setMethod("crs", signature('SpatRaster'), 
	function(x) {
		x@ptr$crs
	}
)


setMethod("crs<-", signature('SpatRaster', 'character'), 
	function(x, ..., value) {
		x@ptr$crs <- value
		return(x)
	}
)


setMethod("crs", signature('SpatLayer'), 
	function(x) {
		x@ptr$crs
	}
)

setMethod("crs<-", signature('SpatLayer', 'character'), 
	function(x, ..., value) {
		x@ptr$crs <- trimws(value)
		return(x)
	}
)
