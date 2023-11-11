# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("names", signature(x="SpatRaster"),
	function(x) {
		nms <- x@cpp$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRaster"),
	function(x, value)  {
		if (is.null(value)) {
			value <- rep("", nlyr(x))
		} else {
			value <- enc2utf8(as.character(value))
			if (length(value) != nlyr(x)) {
				error("names<-", "incorrect number of names")
			}
		}
		x@cpp <- x@cpp$deepcopy()
		if (! x@cpp$setNames(value, FALSE)) {
			error("names<-", "cannot set these names")
		}
		return(x)
	}
)


.names_check <- function(x, value, index, validate, lengthfx) {
	value <- enc2utf8(as.character(value))
	if (!all(index == 1:lengthfx(x))) {
		n <- names(x)
		n[index] <- value
		value <- n
	}
	if (length(value) != lengthfx(x)) {
		error("names<-", "incorrect number of names")
	}
	if (validate) {
		value <- make.names(value, unique = TRUE)
	}
	value
}


setMethod("set.names", signature(x="SpatRaster"),
	function(x, value, index=1:nlyr(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, nlyr)
		if (!x@cpp$setNames(value, FALSE)) {
			error("set.names", "cannot set these names")
		}
		invisible(TRUE)
	}
)


setMethod("names", signature(x="SpatRasterCollection"),
	function(x) {
		nms <- x@cpp$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRasterCollection"),
	function(x, value) {
		x@cpp <- x@cpp$deepcopy()
		x@cpp$names <- enc2utf8(as.character(value))
		x
	}
)

setMethod("set.names", signature(x="SpatRasterCollection"),
	function(x, value, index=1:length(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, length)
		x@cpp$names <- value
		invisible(TRUE)
	}
)



setMethod("names", signature(x="SpatRasterDataset"),
	function(x) {
		nms <- x@cpp$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRasterDataset"),
	function(x, value) {
		x@cpp <- x@cpp$deepcopy()
		x@cpp$names <- enc2utf8(as.character(value))
		x
	}
)

setMethod("set.names", signature(x="SpatRasterDataset"),
	function(x, value, index=1:length(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, length)
		x@cpp$names <- value
		invisible(TRUE)
	}
)


setMethod("varnames", signature(x="SpatRasterDataset"),
	function(x) {
		nms <- x@cpp$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("varnames<-", signature(x="SpatRasterDataset"),
	function(x, value) {
		value <- enc2utf8(as.character(value))
		x@cpp <- x@cpp$deepcopy()
		x@cpp$names <- value
		x
	}
)



setMethod("names", signature(x="SpatVector"),
	function(x) {
		nms <- x@cpp$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("names", signature(x="SpatVectorProxy"),
	function(x) {
		nms <- x@cpp$v$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("names<-", signature(x="SpatVector"),
	function(x, value)  {
		if (length(value) != ncol(x)) {
			error("names<-,SpatVector", "incorrect number of names")
		}
		value <- enc2utf8(as.character(value))
		x@cpp <- x@cpp$deepcopy()
		x@cpp$names <- value
		if (any(names(x) != value)) {
			warn("names<-", "some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)

setMethod("set.names", signature(x="SpatVector"),
	function(x, value, index=1:ncol(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, ncol)
		x@cpp$names <- value
		invisible(TRUE)
	}
)

setMethod("varnames", signature(x="SpatRaster"),
	function(x) {
		nms <- x@cpp$get_sourcenames()
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("varnames<-", signature(x="SpatRaster"),
	function(x, value)  {
		value <- enc2utf8(as.character(value))
		x@cpp <- x@cpp$deepcopy()
		if (!x@cpp$set_sourcenames(value)) {
			error("varnames<-,SpatRaster", "cannot set these names")
		}
		return(x)
	}
)


setMethod("longnames", signature(x="SpatRasterDataset"),
	function(x) {
		nms <- x@cpp$long_names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("longnames", signature(x="SpatRaster"),
	function(x) {
		nms <- x@cpp$get_sourcenames_long()
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("longnames<-", signature(x="SpatRasterDataset"),
	function(x, value)  {
		x@cpp <- x@cpp$deepcopy()
		x@cpp$long_names <- enc2utf8(as.character(value))
		return(x)
	}
)


setMethod("longnames<-", signature(x="SpatRaster"),
	function(x, value)  {
		x@cpp <- x@cpp$deepcopy()
		value <- enc2utf8(as.character(value))
		if (!x@cpp$set_sourcenames_long(value)) {
			error("longnames<-,SpatRaster", "cannot set these names")
		}
		return(x)
	}
)


setMethod("names", signature(x="SpatVectorCollection"),
	function(x) {
		nms <- x@cpp$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatVectorCollection"),
	function(x, value) {
		x@cpp <- x@cpp$deepcopy()
		x@cpp$setNames(enc2utf8(as.character(value)), FALSE)
		x
	}
)
