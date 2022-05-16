
setMethod ("setCats", "SpatRaster", 
	function(x, ...) {
		warn("setCats", "this function will be removed. You can use 'set.cats' instead")
		set.cats(x, ...)
	}
)

## spatstat conflicts

if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}
setMethod("area", signature(x="SpatRaster"), 
	function(x, ...) {
		error("area was removed. Use cellSize or expanse")
	}
)

