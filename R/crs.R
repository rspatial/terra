# Author: Robert J. Hijmans
# Date :  January 2009
# Version 0.9
# License GPL v3


is.proj <- function(crs) {
	substr(crs, 1, 6) == "+proj="
}

.check_proj4_datum <- function(crs) {
	crs <- trimws(tolower(crs))
	if (!is.proj(crs)) return() 
	x <- trimws(unlist(strsplit(crs, "\\+")))
	d <- grep("datum=", x, value=TRUE)
	if (length(d) > 0) {  
		d <- gsub("datum=", "", d)
		if (!(d %in% c("wgs84", "nad83"))) {
			warn("crs", "Only the WGS84 or NAD83 datum can be used in a PROJ.4 string; use WKT2 or authority:code instead")
		}		
	}
	d <- grep("towgs84=", x, value=TRUE)
	if (length(d) > 0) {  
		warn("crs", "+towgs84 parameters in a PROJ4 string are ignored")
	}
}


.proj4 <- function(x) {
	x@ptr$get_crs("proj4")
}

.srs_describe <- function(srs) {
	info <- .SRSinfo(srs)
	names(info) <- c("name", "EPSG", "area", "extent")
	d <- data.frame(t(info), stringsAsFactors=FALSE)
	d$area <- gsub("\\.$", "", d$area)
	d[d == ""] <- NA
	if (is.na(d$extent)) {
		d$extent <- list(c(NA, NA, NA, NA))
	} else {
		d$extent <- list(as.numeric(unlist(strsplit(d$extent, ","))))
	}
	d
}

setMethod("crs", signature("SpatRaster"), 
	function(x, proj=FALSE, describe=FALSE) {
		if (describe) {
			d <- .srs_describe(x@ptr$get_crs("wkt"))
			if (proj) {
				d$proj <- x@ptr$get_crs("proj4")		
			}
			d
		} else if (proj) {
			x@ptr$get_crs("proj4")		
		}  else {
			x@ptr$get_crs("wkt")
		}
	}
)


.txtCRS <- function(x, warn=TRUE) {
	if (inherits(x, "SpatVector") | inherits(x, "SpatRaster")) {
		x <- crs(x)
	}
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
	.check_proj4_datum(x)
	x
}

setMethod("crs<-", signature("SpatRaster", "ANY"), 
	function(x, value) {
		value <- .txtCRS(value)
		x@ptr <- x@ptr$deepcopy()
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
	function(x, proj=FALSE, describe=FALSE) {
		if (describe) {
			d <- .srs_describe(x@ptr$get_crs("wkt"))
			if (proj) {
				d$proj <- x@ptr$get_crs("proj4")		
			}
			d
		} else if (proj) {
			x@ptr$get_crs("proj4")		
		} else {
			x@ptr$get_crs("wkt")
		}
	}
)

setMethod("crs<-", signature("SpatVector", "ANY"), 
	function(x, value) {
		value <- .txtCRS(value)
		x@ptr <- x@ptr$deepcopy()
		x@ptr$set_crs(value)
		messages(x, "crs<-")
	}
)



setMethod("is.lonlat", signature("SpatRaster"), 
	function(x, perhaps=FALSE, warn=TRUE, global=FALSE) {
		if (perhaps) {
			ok <- x@ptr$isLonLat()
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
			ok <- x@ptr$isLonLat()
			if (ok && global) {
				ok <- x@ptr$isGlobalLonLat()
			}
			return(ok)
		}
	}
)


setMethod("is.lonlat", signature("SpatVector"), 
	function(x, perhaps=FALSE, warn=TRUE) {
		ok <- x@ptr$isLonLat()
		if (ok) return(ok)
		if (perhaps) {
			ok <- x@ptr$couldBeLonLat()
			if (ok && warn) {
				if (crs(x) == "") warn("is.lonlat", "assuming lon/lat crs")
			}
		}
		ok
	}
)

