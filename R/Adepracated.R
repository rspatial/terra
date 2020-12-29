desc <- function(x, ...) {
	error("desc", "deprecated function. Use 'describe'")
}

setMethod("isLonLat", signature("SpatRaster"), 
	function(x, ...) {
		warn("isLonLat", "deprecated  method. Use 'is.lonlat'")
		is.lonlat(x, ...)
	}
)
