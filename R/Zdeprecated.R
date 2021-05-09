
## spatstat conflicts

if (!isGeneric("convexhull")) {setGeneric("convexhull", function(x, ...) standardGeneric("convexhull"))}
setMethod("convexhull", signature(x="SpatVector"), 
	function(x) {
		error("convexhull", "terra::convexhull has been removed. Use 'convHull' instead")
	}
)


if (!isGeneric("perimeter")) {setGeneric("perimeter", function(x, ...) standardGeneric("perimeter"))}
setMethod("perimeter", signature(x="SpatVector"), 
	function(x) {
		error("perimeter", "terra::perimeter has been removed. Use 'perim' instead")
	}
)

if (!isGeneric("tiles")) {setGeneric("tiles", function(x,...) standardGeneric("tiles"))}
setMethod("tiles", signature(x="SpatRaster"), 
	function(x) {
		error("tiles", "terra::tiles has been removed. Use 'makeTiles' instead")
	}
)



## tidyverse conflicts
#not exported
if (!isGeneric("expand")) {setGeneric("expand", function(x, y, ...) standardGeneric("expand"))}
setMethod("expand", signature(x="ANY"), 
	function(...) {
		warn("expand", "terra::expand has been removed. Use 'extend' instead")
		extend(...)
	}
)

#not exported
if (!isGeneric("near")) {setGeneric("near", function(x, ...) standardGeneric("near"))}
setMethod("near", signature(x="ANY"), 
	function(x, ...) {
		error("near", "near is deprecated. Use 'nearby' instead")
	}
)

#not exported
if (!isGeneric("separate")) {setGeneric("separate", function(x, ...) standardGeneric("separate"))}
setMethod("separate", signature(x="ANY"), 
	function(x, ...) {
		error("separate", "deprecated function. Use 'segregate'")
		segregate(x, ...)
	}
)


#not exported
if (!isGeneric("pack")) {setGeneric("pack", function(x, ...) standardGeneric("pack"))}
setMethod("pack", signature(x="ANY"), 
	function(x, ...) {
		error("pack", "deprecated function, use 'wrap' instead")
		wrap(x, ...)
	}
)

## replaced for other reasons
#not exported
desc <- function(x) {#
	error("desc", "deprecated function. Use 'describe'")
}


#not exported
gdal_version <- function() {
	gdal()
}

