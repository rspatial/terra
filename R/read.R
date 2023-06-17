# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3

readAll <- function(x) {
	ok <- x@ptr$readAll()
	x <- messages(x, "readAll")
	invisible(ok)
}

setMethod("readStart", signature(x="SpatRaster"),
	function(x) {
		success <- x@ptr$readStart()
		messages(x, "readStart")
		if (!success) error("readStart,SpatRaster", "cannot open file for reading")
		invisible(success)
	}
)

setMethod("readStart", signature(x="SpatRasterDataset"),
	function(x) {
		success <- x@ptr$readStart()
		messages(x, "readStart")
		if (!success) error("readStart,SpatRasterDataset", "cannot open file for reading")
		invisible(success)
	}
)

#setMethod("readStart", signature(x="SpatRasterDataset"),
#	function(x, ...) {
#		nsd <- length(x)
#		for (i in 1:nsd) {
#			y <- x[i]
#			success <- readStart(y)
#			x[i] <- y
#		}
#		messages(x, "readStart")
#		invisible(success)
#	}
#)


setMethod("readStop", signature(x="SpatRaster"),
	function(x) {
		success <- x@ptr$readStop()
		messages(x, "readStop")
		invisible(success)
	}
)

setMethod("readStop", signature(x="SpatRasterDataset"),
	function(x) {
		success <- x@ptr$readStop()
		messages(x, "readStop")
		invisible(success)
	}
)

