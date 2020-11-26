# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("time", signature(x="SpatRaster"), 
	function(x, ...) { 
		#d <- x@ptr$time
		d <- x@ptr$timedbl
		if (x@ptr$timestep == "seconds") {
			strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d
		} else if (x@ptr$timestep == "days") {
			as.Date("1970-01-01") + d
		} else {
			stop("yikes")
		}
	}
)


setMethod("time<-", signature(x="SpatRaster"), 
	function(x, value)  {
		if (length(value) != nlyr(x)) {
			stop("incorrect number of values")
			#if (! x@ptr$setTime(value, "days")) {
			#	stop("cannot set these  values")
			#}
		}
		if (inherits(value, "Date")) {
			value <- as.POSIXlt(value)
		}
		if (inherits(value, "POSIXt")) {
			value <- as.numeric(value)
			if (!x@ptr$setTime(value, "seconds")) {
				stop("cannot set these  values")
			}
		}
		return(x)
	}
)



setMethod("depth", signature(x="SpatRaster"), 
	function(x, ...) { 
		d <- x@ptr$depth
	}
)


setMethod("depth<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- as.numeric(value)
		if (length(value) != nlyr(x)) {
			stop("incorrect number of values")
		}
		if (! x@ptr$setDepth(value)) {
			stop("cannot set these  values")
		}
		#if (any(time(x) != value)) {
		#	warning("some values were changed to make them valid and/or unique")
		#}
		return(x)
	}
)


