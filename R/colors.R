

setMethod ("cols" , "SpatRaster", 
	function(x) {
		hascols <- x@ptr$hasColors()
		if (any(hascols)) {
			d <- x@ptr$getColors()
			d <- lapply(d, .getSpatDF)
		} else {
			d <- vector("list", length(hascols))
		}
		d
	}
)

