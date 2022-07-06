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
		if (!(d %in% c("wgs84", "nad83", "nad27"))) {
			warn("crs<-", "Only the WGS84, NAD83 and NAD27 datums can be used with a PROJ.4 string. Use WKT2, authority:code, or +towgs84= instead")
		}
	}
	#d <- grep("towgs84=", x, value=TRUE)
	#if (length(d) > 0) {  
	#	warn("crs<-", "+towgs84 parameters in a PROJ4 string are ignored")
	#}
}


.proj4 <- function(x) {
	x@ptr$get_crs("proj4")
}

.name_or_proj4 <- function(x) {
	if (inherits(x, "SpatVectorProxy")) {
		v <- vect()
		v@ptr <- x@ptr$v
		x <- v
	}
	d <- .srs_describe(x@ptr$get_crs("wkt"))
	r <- x@ptr$get_crs("proj4")
	if (!(d$name %in% c(NA, "unknown", "unnamed"))) {
		if (substr(r, 1, 13) == "+proj=longlat") {
			r <- paste("lon/lat", d$name)
		} else {
			r <- d$name
		}
		if (!is.na(d$code)) {
			r <- paste0(r, " (", d$authority, ":", d$code, ")")
		} 
	} 
	r
}



.srs_describe <- function(srs) {
	info <- .SRSinfo(srs)
	names(info) <- c("name", "authority", "code", "area", "extent")
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



.get_CRS <- function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
	if (describe) {
		d <- .srs_describe(x@ptr$get_crs("wkt"))
		if (proj) {
			d$proj <- x@ptr$get_crs("proj4")
		}
		d
	} else if (proj) {
		x@ptr$get_crs("proj4")
	} else {
		r <- x@ptr$get_crs("wkt")
		if (parse) {
			unlist(strsplit(r, "\n"))
		} else {
			r
		}
	}
}


setMethod("crs", signature("SpatRaster"), 
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		.get_CRS(x, proj=proj, describe=describe, parse=parse)
	}
)

setMethod("crs", signature("SpatRasterDataset"), 
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		if (length(x) > 0) {
			.get_CRS(x[[1]], proj=proj, describe=describe, parse=parse)
		} else {
			NULL
		}
	}
)


.txtCRS <- function(x, warn=TRUE) {
	if (inherits(x, "SpatVector") | inherits(x, "SpatRaster")) {
		x <- crs(x)
	}
	if (is.null(x) || is.na(x)) {
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

setMethod("set.crs", signature("SpatRaster"), 
	function(x, value) {
		value <- .txtCRS(value)
		x@ptr$set_crs(value)
		messages(x, "set_crs")
	}
)



#setMethod("crs<-", signature("SpatRaster", "character"), 
#	function(x, ..., value) {
#		x@ptr$set_crs(value[1])
#		messages(x, "crs<-")
#	}
#)


setMethod("crs", signature("SpatVector"), 
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		.get_CRS(x, proj=proj, describe=describe, parse=parse)
	}
)

setMethod("crs", signature("SpatVectorProxy"), 
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		v <- vect()
		v@ptr <- x@ptr$v
		.get_CRS(v, proj=proj, describe=describe, parse=parse)
	}
)

setMethod("crs", signature("sf"),
  function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
	  xcrs <- attr(x[[ attr(x, "sf_column") ]], "crs")$wkt
	  x <- vect(cbind(0,0), crs=xcrs)
	  .get_CRS(x, proj=proj, describe=describe, parse=parse)
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

setMethod("set.crs", signature("SpatVector"), 
	function(x, value) {
		value <- .txtCRS(value)
		x@ptr$set_crs(value)
		messages(x, "set_crs")
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
		} else {
			ok <- x@ptr$isLonLat()
			if (ok && global) {
				ok <- x@ptr$isGlobalLonLat()
			}
			if ((!ok) && (crs(x) == "")) {
				ok <- NA
			}
		}
		ok
	}
)


setMethod("is.lonlat", signature("SpatVector"), 
	function(x, perhaps=FALSE, warn=TRUE) {
		if (perhaps) {
			ok <- x@ptr$couldBeLonLat()
			if (ok && warn) {
				if (crs(x) == "") warn("is.lonlat", "assuming lon/lat crs")
			}
		} else {
			ok <- x@ptr$isLonLat()
			if ((!ok) && (crs(x) == "")) {
				ok <- NA
			}
		}
		ok
	}
)

