# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("time", signature(x="SpatRaster"), 
	function(x, ...) { 
		if (!x@ptr$hasTime) {
			return(rep(NA, nlyr(x)))
		}
		d <- x@ptr$time
		tstep <- x@ptr$timestep 
		if (tstep == "seconds") {
			strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d
		} else if (tstep == "days") {
			as.Date("1970-01-01") + d
		#} else if (tstep == "months") {
		#} else if (tstep == "years") {
		} else { # raw 
			d
		}
	}
)


setMethod("time<-", signature(x="SpatRaster"), 
	function(x, value)  {
		if (inherits(value, "Date")) {
			value <- as.POSIXlt(value)
		}
		if (inherits(value, "POSIXt")) {
			if (!x@ptr$setTime(as.numeric(value), "seconds")) {
				error("time<-", "cannot set these  values")
			}
		} else {
			if (!x@ptr$setTime(as.numeric(value), "raw")) {
				error("time<-", "cannot set these  values")
			}
		}
		return(x)
	}
)



setMethod("depth", signature(x="SpatRaster"), 
	function(x, ...) { 
		x@ptr$depth
	}
)


setMethod("depth<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- as.numeric(value)
		if (! x@ptr$setDepth(value)) {
			error("depth<-", "cannot set these  values")
		}
		return(x)
	}
)

setMethod("linearUnits", signature(x="SpatRaster"), 
	function(x, ...) {
		.getLinearUnits(crs(x))
	}
)

setMethod("linearUnits", signature(x="SpatVector"), 
	function(x, ...) {
		.getLinearUnits(crs(x))
	}
)

setMethod("units", signature(x="SpatRaster"), 
	function(x) { 
		x@ptr$units
	}
)

setMethod("units<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- as.character(value)
		if (! x@ptr$set_units(value)) {
			error("units<-", "cannot set these  values")
		}
		return(x)
	}
)


setMethod("units", signature(x="SpatRasterDataset"), 
	function(x) { 
		x@ptr$units
	}
)

setMethod("units<-", signature(x="SpatRasterDataset"), 
	function(x, value)  {
		value <- as.character(value)
		x@ptr$units <- value
		return(x)
	}
)

