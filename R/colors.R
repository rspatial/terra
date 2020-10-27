

setMethod ("cols" , "SpatRaster", 
	function(x) {
		hascols <- x@ptr$hasColors()
		if (any(hascols)) {
			d <- x@ptr$getColors()
			d <- lapply(d, terra:::.getSpatDF)
		} else {
			d <- vector("list", length(att))
		}
		d
	}
)

