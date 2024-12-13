# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("names", signature(x="SpatRaster"),
	function(x) {
		nms <- x@pntr$names
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
		x@pntr <- x@pntr$deepcopy()
		if (! x@pntr$setNames(value, FALSE)) {
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
		if (!x@pntr$setNames(value, FALSE)) {
			error("set.names", "cannot set these names")
		}
		invisible(TRUE)
	}
)


setMethod("names", signature(x="SpatRasterCollection"),
	function(x) {
		nms <- x@pntr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRasterCollection"),
	function(x, value) {
		x@pntr <- x@pntr$deepcopy()
		if (is.null(value)) {
			value <- rep("", length(x))
		}
		x@pntr$names <- enc2utf8(as.character(value))
		x
	}
)

setMethod("set.names", signature(x="SpatRasterCollection"),
	function(x, value, index=1:length(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, length)
		x@pntr$names <- value
		invisible(TRUE)
	}
)



setMethod("names", signature(x="SpatRasterDataset"),
	function(x) {
		nms <- x@pntr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRasterDataset"),
	function(x, value) {
		x@pntr <- x@pntr$deepcopy()
		if (is.null(value)) {
			value <- rep("", length(x))
		}
		if (is.list(value)) {
			nl <- nlyr(x)
			if (length(value) == 1) {
				if (length(unique(nl)) > 1) {
					error("names<-", "the number of layers varies between datasets")
				}
				x@pntr$set_layernames(enc2utf8(as.character(value[[1]])), -1)
			} else {
				if (length(value) != length(x)) {
					error("names<-", "the number of list elements does not match the number of datasets")				
				}
				for (i in seq_along(length(x))) x@pntr$set_layernames(enc2utf8(as.character(value[[i]])), i-1)
			}
		} else {
			x@pntr$names <- enc2utf8(as.character(value))
		}
		messages(x, "names<-")
	}
)


setMethod("set.names", signature(x="SpatRasterDataset"),
	function(x, value, index=1:length(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, length)
		x@pntr$names <- value
		invisible(TRUE)
	}
)


setMethod("varnames", signature(x="SpatRasterDataset"),
	function(x) {
		nms <- x@pntr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("varnames<-", signature(x="SpatRasterDataset"),
	function(x, value) {
		value <- enc2utf8(as.character(value))
		x@pntr <- x@pntr$deepcopy()
		x@pntr$names <- value
		x
	}
)



setMethod("names", signature(x="SpatVector"),
	function(x) {
		nms <- x@pntr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("names", signature(x="SpatVectorProxy"),
	function(x) {
		nms <- x@pntr$v$names
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
		x@pntr <- x@pntr$deepcopy()
		x@pntr$names <- value
		if (any(names(x) != value)) {
			warn("names<-", "some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)

setMethod("set.names", signature(x="SpatVector"),
	function(x, value, index=1:ncol(x), validate=FALSE)  {
		value <- .names_check(x, value, index, validate, ncol)
		x@pntr$names <- value
		invisible(TRUE)
	}
)

setMethod("varnames", signature(x="SpatRaster"),
	function(x) {
		nms <- x@pntr$get_sourcenames()
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("varnames<-", signature(x="SpatRaster"),
	function(x, value)  {
		value <- enc2utf8(as.character(value))
		x@pntr <- x@pntr$deepcopy()
		if (!x@pntr$set_sourcenames(value)) {
			error("varnames<-,SpatRaster", "cannot set these names")
		}
		return(x)
	}
)


setMethod("longnames", signature(x="SpatRasterDataset"),
	function(x) {
		nms <- x@pntr$long_names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("longnames", signature(x="SpatRaster"),
	function(x) {
		nms <- x@pntr$get_sourcenames_long()
		Encoding(nms) <- "UTF-8"
		nms
	}
)

setMethod("longnames<-", signature(x="SpatRasterDataset"),
	function(x, value)  {
		x@pntr <- x@pntr$deepcopy()
		x@pntr$long_names <- enc2utf8(as.character(value))
		return(x)
	}
)


setMethod("longnames<-", signature(x="SpatRaster"),
	function(x, value)  {
		x@pntr <- x@pntr$deepcopy()
		value <- enc2utf8(as.character(value))
		if (!x@pntr$set_sourcenames_long(value)) {
			error("longnames<-,SpatRaster", "cannot set these names")
		}
		return(x)
	}
)


setMethod("names", signature(x="SpatVectorCollection"),
	function(x) {
		nms <- x@pntr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatVectorCollection"),
	function(x, value) {
		x@pntr <- x@pntr$deepcopy()
		if (is.null(value)) {
			value <- rep("", length(x))
		}
		x@pntr$setNames(enc2utf8(as.character(value)), FALSE)
		x
	}
)
