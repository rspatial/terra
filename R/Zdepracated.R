
setMethod("near", signature(x="ANY"), 
	function(x, ...) {
		warn("near", "near is deprecated. Use 'nearby' instead")
		nearby(x, ...)
	}
)


setMethod("separate", signature(x="ANY"), 
	function(x, ...) {
		warn("separate", "deprecated function. Use 'segregate'")
		segregate(x, ...)
	}
)


setMethod("pack", signature(x="ANY"), 
	function(x, ...) {
		warn("pack", "deprecated function, use 'wrap' instead")
		wrap(x, ...)
	}
)


desc <- function(x) {
	error("desc", "deprecated function. Use 'describe'")
}



gdal_version <- function() {
	gdal()
}


