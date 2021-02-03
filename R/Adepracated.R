desc <- function(x) {
	error("desc", "deprecated function. Use 'describe'")
}

setMethod("isLonLat", signature("SpatRaster"), 
	function(x) {
		error("isLonLat", "deprecated  method. Use 'is.lonlat'")
		#is.lonlat(x)
	}
)
