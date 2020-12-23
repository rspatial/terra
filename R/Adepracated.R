desc <- function(x, ...) {
	error("desc", "depracated function. Use 'describe'")
}

setMethod("isLonLat", signature("SpatRaster"), 
	function(x, ...) {
		warn("isLonLat", "depracated  method. Use 'is.lonlat'")
		is.lonlat(x, ...)
	}
)
