# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("time", signature(x="SpatRaster"), 
	function(x, format="") { 
		if (!x@ptr$hasTime) {
			return(rep(NA, nlyr(x)))
		}
		d <- x@ptr$time
		if (format != "") {
			tstep <- format
		} else {
			tstep <- x@ptr$timestep 
		}
		if (tstep == "seconds") {
			d <- strptime("1970-01-01", "%Y-%m-%d", tz="UTC") + d
			tz <- x@ptr$timezone
			if (!(tz %in% c("", "UTC"))) {
				attr(d, "tzone") = tz
			}
			d
		} else if (tstep == "days") {
			d <- strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d
			as.Date(d)
		} else if (tstep == "yearmonths") {
			d <- strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d
			y <- as.integer(format(d, "%Y"))
			y + (as.integer(format(d, "%m"))-1)/12
		} else if (tstep == "months") {
			d <- strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d
			as.integer(format(d, "%m"))
		} else if (tstep == "years") {
			d <- strptime("1970-01-01", "%Y-%m-%d", tz = "UTC") + d
			as.integer(format(d, "%Y"))
		} else { # raw 
			d
		}
	}
)


setMethod("time<-", signature(x="SpatRaster"), 
	function(x, value, tstep="")  {
		if (is.null(value)) {
			x@ptr$setTime(0[0], "remove")
			return(x)
		} 
		if (inherits(value, "character")) {
			error("time<-", "value cannot be a character type")
		}
		tzone <- "UTC"
		if (inherits(value, "Date")) {
			value <- as.POSIXlt(value)
			tstep <- "days"
		} else if (inherits(value, "POSIXt")) {
			tstep <- "seconds"
			tzone <- attr(value, "tzone")
			if (is.null(tzone)) tzone = ""
		} else if (inherits(value, "yearmon")) {
			value <- as.numeric(value)
			year <- floor(value)
			month <- round(12 * (value - year) + 1)
			d <- as.Date(paste(year, month, "15", sep="-")) 
			value <- as.POSIXlt(d)
			tstep <- "yearmonths"
		} else if (tstep == "years") {
			value <- as.POSIXlt(as.Date(paste0(floor(value), "-6-15")))
		} else if (tstep == "months") {
			value <- floor(value)
			if (!all(value %in% 1:12)) {
				error("date<-", "month values should be between 1 and 12")
			}
			value <- as.POSIXlt(as.Date(paste0("1970-", value, "-15")))
		} else if (tstep == "") {
			tstep <- "raw"
		} else {
			error("time<-", "unknown tstep")
		}
		if (!x@ptr$setTime(as.numeric(value), tstep, tzone)) {
			error("time<-", "cannot set these values")
		}
		return(x)
	}
)



setMethod("depth", signature(x="SpatRaster"), 
	function(x) { 
		x@ptr$depth
	}
)


setMethod("depth<-", signature(x="SpatRaster"), 
	function(x, value)  {
		if (is.null(value)) {
			x@ptr$setTime(0[0])
			return(x)
		}
		value <- as.numeric(value)
		if (! x@ptr$setDepth(value)) {
			error("depth<-", "cannot set these  values")
		}
		return(x)
	}
)

setMethod("linearUnits", signature(x="SpatRaster"), 
	function(x) {
		.getLinearUnits(crs(x))
	}
)

setMethod("linearUnits", signature(x="SpatVector"), 
	function(x) {
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
		if (is.null(value) || all(is.na(value))) {
			value <- ""
		} else {
			value <- as.character(value)
		}
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

