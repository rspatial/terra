
#not exported
setMethod("expand", signature(x="ANY"), 
	function(...) {
		warn("expand", "terra::expand has been removed. Use 'extend' instead")
		extend(...)
	}
)

#not exported
setMethod("near", signature(x="ANY"), 
	function(x, ...) {
		error("near", "near is deprecated. Use 'nearby' instead")
	}
)

#not exported
setMethod("separate", signature(x="ANY"), 
	function(x, ...) {
		error("separate", "deprecated function. Use 'segregate'")
		segregate(x, ...)
	}
)


#not exported
setMethod("pack", signature(x="ANY"), 
	function(x, ...) {
		error("pack", "deprecated function, use 'wrap' instead")
		wrap(x, ...)
	}
)


#not exported
desc <- function(x) {#
	error("desc", "deprecated function. Use 'describe'")
}


#not exported
gdal_version <- function() {
	gdal()
}


