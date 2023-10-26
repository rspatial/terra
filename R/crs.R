# Author: Robert J. Hijmans
# Date :  January 2009
# Version 0.9
# License GPL v3



character_crs <- function(x, caller="") {

	if (!inherits(x, "character")) {
		if (is.atomic(x)) {
			# for logical NA 
			x <- as.character(x)
		} else {
			x <- crs(x)
		}
	}

	if (is.na(x)) {
		""
	} else {
		if (tolower(x) == "local") {
			x <- 'LOCAL_CS["Cartesian (Meter)", LOCAL_DATUM["Local Datum",0], UNIT["Meter",1.0], AXIS["X",EAST], AXIS["Y",NORTH]]'
		} else if (tolower(x) == "lonlat") {
			x <- "+proj=longlat"
		}
		x
	}
}



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
	x@cpp$get_crs("proj4")
}

.name_from_wkt <- function(wkt) {
	s = strsplit(wkt, ",")[[1]][1]
	strsplit(s, "\"")[[1]][[2]]
}

.name_or_proj4 <- function(x) {
	if (inherits(x, "SpatVectorProxy")) {
		ptr <- x@cpp$v
	} else if (inherits(x, "Rcpp_SpatRaster")) {
		ptr <- x	
	} else {
		ptr <- x@cpp
	}
	wkt <- ptr$get_crs("wkt")
	d <- .srs_describe(wkt)
	r <- ptr$get_crs("proj4")
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
	if (r == "") {
		rr <- try(.name_from_wkt(wkt), silent=TRUE)
		if (!inherits(rr, "try-error")) {
			r <- rr
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
		d <- .srs_describe(x@cpp$get_crs("wkt"))
		if (proj) {
			d$proj <- x@cpp$get_crs("proj4")
		}
		d
	} else if (proj) {
		x@cpp$get_crs("proj4")
	} else {
		r <- x@cpp$get_crs("wkt")
		if (parse) {
			unlist(strsplit(r, "\n"))
		} else {
			r
		}
	}
}


setMethod("crs", signature("character"),
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		x <- rast(crs=x)
		.get_CRS(x, proj=proj, describe=describe, parse=parse)
	}
)

setMethod("crs", signature("SpatExtent"),
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		return("")
	}
)

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
		lowx <- tolower(x)
		if (lowx == "local") {
			x = 'LOCAL_CS["Cartesian (Meter)", LOCAL_DATUM["Local Datum",0], UNIT["Meter",1.0], AXIS["X",EAST], AXIS["Y",NORTH]]'
		} else if (lowx == "lonlat") {
			x <- "+proj=longlat"
		}
	} else {
		error("crs", "I do not know what to do with this argument (expected a character string)")
	}
	.check_proj4_datum(x)
	x
}

setMethod("crs<-", signature("SpatRaster", "ANY"),
	function(x, warn=FALSE, value) {
		if (missing(value)) {
			value <- warn
			warn <- FALSE
		}
		value <- .txtCRS(value)
		if (warn && (crs(x) != "") && (value != "")) {
			message("Assigning a new crs. Use 'project' to transform a SpatRaster to a new crs")
		}
		x@cpp <- x@cpp$deepcopy()
		x@cpp$set_crs(value)
		messages(x, "crs<-")
	}
)

setMethod("set.crs", signature("SpatRaster"),
	function(x, value) {
		value <- .txtCRS(value)
		x@cpp$set_crs(value)
		messages(x, "set_crs")
		invisible(TRUE)
	}
)



#setMethod("crs<-", signature("SpatRaster", "character"),
#	function(x, ..., value) {
#		x@cpp$set_crs(value[1])
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
		v@cpp <- x@cpp$v
		.get_CRS(v, proj=proj, describe=describe, parse=parse)
	}
)

setMethod("crs", signature("SpatVectorCollection"),
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		if (length(x) > 0) {
			.get_CRS(x[[1]], proj=proj, describe=describe, parse=parse)
		} else {
			NULL
		}
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
	function(x, warn=FALSE, value) {
		if (missing(value)) {
			value <- warn
			warn <- FALSE
		}
		value <- .txtCRS(value)
		if (warn && (crs(x) != "") && (value != "")) {
			message("Assigning a new crs. Use 'project' to transform a SpatVector to a new crs")
		}
		x@cpp <- x@cpp$deepcopy()
		x@cpp$set_crs(value)
		messages(x, "crs<-")
	}
)


setMethod("set.crs", signature("SpatVector"),
	function(x, value) {
		value <- .txtCRS(value)
		x@cpp$set_crs(value)
		messages(x, "set_crs")
	}
)


setMethod("is.lonlat", signature("SpatRaster"),
	function(x, perhaps=FALSE, warn=TRUE, global=FALSE) {
		if (perhaps) {
			ok <- x@cpp$isLonLat()
			if (ok) {
				messages(x, "is.lonlat")
				if (global) {
					return(x@cpp$isGlobalLonLat())
				} else {
					return(ok)
				}
			} 
			ok <- x@cpp$couldBeLonLat()
			if (ok) {
				messages(x, "is.lonlat")
				if (global) {
					ok <- x@cpp$isGlobalLonLat()
				}
			}
			if (ok && warn) {
				warn("is.lonlat", "assuming lon/lat crs")
			}
		} else {
			ok <- x@cpp$isLonLat()
			if (ok) {
				messages(x, "is.lonlat")
				if (global) {
					ok <- x@cpp$isGlobalLonLat()
				}
			} else if (crs(x) == "") {
				ok <- NA
				warn("is.lonlat", "unknown crs")
			}
		}
		ok
	}
)


setMethod("is.lonlat", signature("SpatVector"),
	function(x, perhaps=FALSE, warn=TRUE) {
		if (perhaps) {
			ok <- x@cpp$isLonLat()
			if (ok) {
				messages(x, "is.lonlat")
				return(ok)
			} 
			ok <- x@cpp$couldBeLonLat()
			if (ok) {
				messages(x, "is.lonlat")
			}
			if (ok && warn) {
				warn("is.lonlat", "assuming lon/lat crs")
			}
		} else {
			ok <- x@cpp$isLonLat()
			if (ok) {
				messages(x, "is.lonlat")
			} else if (crs(x) == "") {
				ok <- NA
				warn("is.lonlat", "unknown crs")
			}
		}
		ok
	}
)

setMethod("is.lonlat", signature("character"),
	function(x, perhaps=FALSE, warn=TRUE) {
		x <- rast(crs=x)
		is.lonlat(x, perhaps=perhaps, warn=warn)
	}
)

same.crs <- function(x, y) {
	if (!is.character(x)) {
		x <- crs(x)
	}
	if (!is.character(y)) {
		y <- crs(y)
	}
	if (inherits(x, "CRS")) {
           if (!is.null(comment(x))) {
			x <- comment(x)
		} else {
			x <- x@projargs
           } 
	}
	if (inherits(y, "CRS")) {
           if (!is.null(comment(y))) {
			y <- comment(y)
		} else {
			y <- y@projargs
           } 
	}
	if (is.na(x)) x <- ""
	if (is.na(y)) x <- ""
	if (!is.character(x)) {
		x <- as.character(x)
	}
	if (!is.character(y)) {
		y <- as.character(y)
	}
	.sameSRS(x, y)
}



