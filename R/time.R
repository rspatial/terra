# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


format_ym <- function(x) {
	y <- floor(x)
	m <- round((x-y) * 12 + 1)
	m <- month.abb[m]
	paste(y, m, sep="-")
}

yearweek <- function(d) {
	y <- as.integer(strftime(d, format = "%Y"))
	w <- strftime(d, format = "%V")
	m <- strftime(d, format = "%m")
	i <- w > "51" & m=="01"
	y[i] <- y[i] - 1
	i <- w=="01" & m=="12"
	y[i] <- y[i] + 1
	yy <- as.character(y)
	i <- nchar(yy) < 4
	yy[i] <- formatC(y[i], width=4, flag="0")	
	paste0(yy, w)
}

setMethod("has.time", signature(x="SpatRaster"),
	function(x) {
		x@pntr$hasTime
	}
)


setMethod("timeInfo", signature(x="SpatRaster"),
	function(x) {
		time <- x@pntr$hasTime
		if (time) {
			step <- x@pntr$timestep
			if (step == "seconds") {
				data.frame(time=time, step=step, zone=x@pntr$timezone)
			} else {
				data.frame(time=time, step=step, zone="")
			}
		} else {
			data.frame(time=time, step="", zone="")
		}
	}
)


setMethod("timeInfo", signature(x="SpatRasterDataset"),
	function(x) {
		t(sapply(x, timeInfo))
	}
)


time_as_seconds <- function(x) {
	d <- x@pntr$time
	d <- strptime("1970-01-01", "%Y-%m-%d", tz="UTC") + d
	tz <- x@pntr$timezone
	if (!(tz %in% c("", "UTC"))) {
		attr(d, "tzone") = tz
	}
	d
}


#setMethod("time", signature(x="SpatVector"),
#	function(x, format="") {
#		cls <- sapply(values(x[1,]), function(i) { a = class(i); a[length(a)] })
#		i <- which(cls %in% c("Date", "POSIXt"))[1]
#		if (is.na(i)) {
#			return(rep(NA, nrow(x)))
#		} else {
#			d <- x[,i,drop=TRUE][,,drop=TRUE]
#			if (format != "") {
#				steps <- c("seconds", "days", "months", "years", "yearmonths")
#				format <- match.arg(tolower(format), steps)
#				if (!(format %in% steps)) {
#					error("time", "not a valid time format")
#				}
#				tstep <- ifelse(cls[i]=="Date", "days", "seconds")
#				if (format == "seconds") {
#					if (tstep != "seconds") {
#						error("time", "cannot extract seconds from Dates")
#					}
#					d
#				} else if (format == "days") {
#					as.Date(d)
#				} else if (format == "yearmonths") {
#					y <- as.integer(format(d, "%Y"))
#					y + (as.integer(format(d, "%m"))-1)/12
#				} else if (format == "months") {
#					as.integer(format(d, "%m"))
#				} else if (format == "years") {
#					as.integer(format(d, "%Y"))
#				}
#			} else {
#				d
#			}
#		}
#	}
#)

setMethod("time", signature(x="SpatRaster"),
	function(x, format="") {
		if (!x@pntr$hasTime) {
			return(rep(NA, nlyr(x)))
		}
		d <- x@pntr$time
		tstep <- x@pntr$timestep
		
		if (format != "") {
			steps <- c("seconds", "days", "months", "years", "yearmonths")
			format <- match.arg(tolower(format), steps)
			if ((format == "months") && (tstep == "years")) {
				error("time", "cannot extract months from years-time")
			} else if ((format == "years") && (tstep %in% c("months"))) {
				error("time", "cannot extract years from months-time")
			} else if ((format == "yearmonths") && (tstep %in% c("months", "years"))) {
				error("time", "cannot extract yearmonths from this type of time data")
			} else if ((format == "seconds") && (tstep != "seconds")) {
				error("time", "cannot extract seconds from this type of time data")
			} else if ((format == "days") && (!(tstep %in% c("seconds", "days")))) {
				error("time", "cannot extract days from this type of time data")
			}
			tstep <- format
		} else if (tstep == "raw") {
			return(d)
		}
		
		
		d <- strptime("1970-01-01", "%Y-%m-%d", tz="UTC") + d
		if (tstep == "seconds") {
			tz <- x@pntr$timezone
			if (!(tz %in% c("", "UTC"))) {
				attr(d, "tzone") = tz
			}
			d
		} else if (tstep == "days") {
			as.Date(d)
		} else if (tstep == "yearmonths") {
			y <- as.integer(format(d, "%Y"))
			y + (as.integer(format(d, "%m"))-1)/12
		} else if (tstep == "months") {
			as.integer(format(d, "%m"))
		} else if (tstep == "years") {
			as.integer(format(d, "%Y"))
#		} else if (tstep == "yearweeks") {
#			yearweek(as.Date(d))
		} else { # ???
			d
		}
	}
)

setMethod("time", signature(x="SpatRasterDataset"),
	function(x, format="") {
		lapply(x,time, format=format)
	}
)

posix_from_ym <- function(y, m) {
	y <- floor(y)
	i <- ((y < 0) | (y > 9999))
	if (any(i)) {
		d <- paste(paste(rep("1900", length(y)), m, "15", sep="-"), "12:00:00")
		d[!i] <- paste(paste(y[!i], m, "15", sep="-"), "12:00:00")
		d <- as.POSIXlt(d, format="%Y-%m-%d %H:%M:%S", tz="UTC")
		for (j in i) {
			d$year[j] = y[j] - 1900
		}
		d
	} else {
		d <- paste(paste(y, m, "15", sep="-"), "12:00:00")
		as.POSIXlt(d, format="%Y-%m-%d %H:%M:%S", tz="UTC")
	}
}


setMethod("time<-", signature(x="SpatRaster"),
	function(x, tstep="", value)  {
		if (missing(value)) {
			value <- tstep
			tstep <- ""
		}
		x@pntr <- x@pntr$deepcopy()
		if (is.null(value)) {
			x@pntr$setTime(0[0], "remove", "")
			return(x)
		}
		if (inherits(value, "character")) {
			error("time<-", "value cannot be a character type")
		}
		if (length(value) != nlyr(x)) {
			error("time<-", "length(value) != nlyr(x)")
		}
		if (tstep != "") {
			tstep = match.arg(as.character(tstep), c("seconds", "days", "months", "years", "yearmonths", "raw"))
		}
		## may not be necessary
		if (tstep == "seconds") tstep = ""
		
		tzone <- "UTC"
		stept <- ""
		if (inherits(value, "Date")) {
			value <- as.POSIXlt(value)
			if (tstep == "") stept <- "days"
		} else if (inherits(value, "POSIXt")) {
			if (tstep == "") stept <- "seconds"
			tzone <- attr(value, "tzone")[1]
			if (is.null(tzone)) tzone = "UTC"
		} else if (inherits(value, "yearmon")) {
			value <- as.numeric(value)
			year <- floor(value)
			month <- round(12 * (value - year) + 1)
			value <- posix_from_ym(value, month)
			if (tstep == "") stept <- "yearmonths"
		} 
		
		if (stept == "") {
			stept = tstep
			if (tstep == "years") {
				if (is.numeric(value)) {
					value <- posix_from_ym(value, "6")
				} else {
					value <- as.integer(strftime(value, format = "%Y", tz=tzone))
					value <- posix_from_ym(value, "6")
				}
			} else if (tstep == "months") {
				if (is.numeric(value)) {
					value <- floor(value)
				} else {
					value <- as.integer(strftime(value, format = "%m", tz=tzone))
				}
				if (!all(value %in% 1:12)) {
					error("date<-", "months should be between 1 and 12")
				}
				value <- posix_from_ym(1970, value)
			} else if (tstep == "yearmonths") {
				if (is.numeric(value)) {
					y <- as.integer(substr(value, 1, 4))
					m <- value - (y * 100)
				} else {
					y <- as.integer(strftime(value, format = "%Y", tz=tzone))
					m <- as.integer(strftime(value, format = "%m", tz=tzone))
				}
				if (!all(m %in% 1:12)) {
					error("date<-", "months should be between 1 and 12")
				}
				value <- posix_from_ym(y, m)
			#} else if (tstep == "days") {
			#	print(value)
			#	value <- as.Date(value)
			#	stept = tstep
			} else if (tstep == "") {
				stept <- "raw"
			}
		}
		if (!x@pntr$setTime(as.numeric(value), stept, tzone)) {
			error("time<-", "cannot set these values")
		}
		return(x)
	}
)


setMethod("time<-", signature(x="SpatRasterDataset"),
	function(x, tstep="", value)  {

		if (missing(value)) {
			value <- tstep
			tstep <- ""
		}
		tstep <- rep_len(tstep, length(x))

		if (is.list(value)) {
			if (length(x) != length(value)) {
				error("time<-", "the list should have the same length as 'x'")
			}
			z <- lapply(1:length(x), function(i) { 
				time(x[i], tstep=tstep[i]) <- value[[i]]
			})
			
		} else {
			if (length(unique(nlyr(x))) > 1) {
				error("time<-", "not all SpatRasters have the same number of layers")
			}
			z <- lapply(1:length(x), function(i) { 
				time(x[i], tstep=tstep[i]) <- value
			})
		}
		x
	}
)


setMethod("depth", signature(x="SpatRaster"),
	function(x) {
		x@pntr$depth
	}
)


setMethod("depth<-", signature(x="SpatRaster"),
	function(x, value)  {
		if (is.null(value)) {
			x@pntr$setDepth(0[0])
			return(x)
		}
		value <- as.numeric(value)
		if (! x@pntr$setDepth(value)) {
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
		x@pntr$units
	}
)

setMethod("units<-", signature(x="SpatRaster"),
	function(x, value)  {
		if (is.null(value) || all(is.na(value))) {
			value <- ""
		} else {
			value <- as.character(value)
		}
		if (! x@pntr$set_units(value)) {
			error("units<-", "cannot set these  values")
		}
		return(x)
	}
)


setMethod("units", signature(x="SpatRasterDataset"),
	function(x) {
		x@pntr$units
	}
)

setMethod("units<-", signature(x="SpatRasterDataset"),
	function(x, value)  {
		value <- as.character(value)
		x@pntr$units <- value
		return(x)
	}
)

