# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


setMethod("names", signature(x="SpatRaster"), 
	function(x) { 
		nms <- x@ptr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- enc2utf8(as.character(value))
		if (length(value) != nlyr(x)) {
			error("names<-", "incorrect number of names")
		}
		if (! x@ptr$setNames(value, FALSE)) {
			error("names<-", "cannot set these names")
		}

		if (any(names(x) != value)) {
			# should only be possible with $setNames(value, TRUE)
			warn("names<-", "some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)

setMethod("names", signature(x="SpatRasterDataset"), 
	function(x) { 
		nms <- x@ptr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("names<-", signature(x="SpatRasterDataset"), 
	function(x, value) {
		x@ptr$names <- enc2utf8(as.character(value))
		x
	}
)

setMethod("varnames", signature(x="SpatRasterDataset"), 
	function(x) {
		nms <- x@ptr$names
		Encoding(nms) <- "UTF-8"
		nms
	}
)


setMethod("varnames<-", signature(x="SpatRasterDataset"), 
	function(x, value) {
		value <- enc2utf8(as.character(value))
		x@ptr$names <- value
		x
	}
)


setMethod("names", signature(x="SpatVector"), 
	function(x) { 
		nms <- x@ptr$names
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
		x@ptr$names <- value
		if (any(names(x) != value)) {
			warn("names<-", "some names were changed to make them valid and/or unique")
		}
		return(x)
	}
)


setMethod("varnames", signature(x="SpatRaster"), 
	function(x) { 
		nms <- x@ptr$get_sourcenames()
		Encoding(nms) <- "UTF-8"
		nms		
	}
)

setMethod("varnames<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- enc2utf8(as.character(value))
		if (!x@ptr$set_sourcenames(value)) {
			error("varnames<-,SpatRaster", "cannot set these names")
		}
		return(x)
	}
)

setMethod("longnames", signature(x="SpatRasterDataset"), 
	function(x) { 
		nms <- x@ptr$long_names
		Encoding(nms) <- "UTF-8"
		nms		
	}
)


setMethod("longnames", signature(x="SpatRaster"), 
	function(x) { 
		nms <- x@ptr$get_sourcenames_long()
		Encoding(nms) <- "UTF-8"
		nms		
	}
)

setMethod("longnames<-", signature(x="SpatRasterDataset"), 
	function(x, value)  {
		x@ptr$long_names <- enc2utf8(as.character(value))
		return(x)
	}
)


setMethod("longnames<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- enc2utf8(as.character(value))
		if (!x@ptr$set_sourcenames_long(value)) {
			error("longnames<-,SpatRaster", "cannot set these names")
		}
		return(x)
	}
)

