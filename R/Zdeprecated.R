
#if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}
#setMethod ("area" , "SpatRaster",
#	function (x, ...) {
#		error("area", "this method was removed. Use expanse or cellSize")
#	}
#)


if (!isGeneric("setCats")) { setGeneric("setCats", function(x, ...) standardGeneric("setCats")) }

setMethod ("setCats" , "SpatRaster",
	function (x, ...) {
		warn("setCats", "this function will be removed. Please can use 'levels<-' or 'set.cats' instead")
		set.cats(x, ...)
	}
)
