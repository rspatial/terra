
## spatstat conflicts

if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}

setMethod("area", signature(x="SpatRaster"), 
	function(x, sum=TRUE, correct=FALSE, mask=FALSE, filename="", ...) {
		if (!sum) {
			error("area", 'area(x, sum=FALSE) will be removed. Use "cellSize(x)"')
			#cellSize(x, mask=mask, transform=correct, filename=filename, ...)
		} else {
			error("area", 'area(x, sum=TRUE) will be removed. Use "expanse(x)" or "global(cellSize(x), "sum")"')
			#expanse(x, transform=correct, ...)			
		}		
	}
)

setMethod("area", signature(x="SpatVector"), 
	function(x, ...) {
		error("area", 'area was removed. Use "expanse(x)"')
		#expanse(x)
	}
)


## old not exported
################
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
if (!isGeneric("expand")) {setGeneric("expand", function(x, y, ...) standardGeneric("expand"))}
setMethod("expand", signature(x="ANY", y="ANY"), 
	function(...) {
		warn("expand", "terra::expand has been removed. Use 'extend' instead")
		extend(...)
	}
)

if (!isGeneric("near")) {setGeneric("near", function(x, ...) standardGeneric("near"))}
setMethod("near", signature(x="ANY"), 
	function(x, ...) {
		error("near", "near is deprecated. Use 'nearby' instead")
	}
)

if (!isGeneric("separate")) {setGeneric("separate", function(x, ...) standardGeneric("separate"))}
setMethod("separate", signature(x="ANY"), 
	function(x, ...) {
		error("separate", "deprecated function. Use 'segregate'")
		segregate(x, ...)
	}
)


if (!isGeneric("pack")) {setGeneric("pack", function(x, ...) standardGeneric("pack"))}
setMethod("pack", signature(x="ANY"), 
	function(x, ...) {
		error("pack", "deprecated function, use 'wrap' instead")
		wrap(x, ...)
	}
)
