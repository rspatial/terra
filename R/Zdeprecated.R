
setMethod ("setCats", "SpatRaster", 
	function(x, ...) {
		warn("setCats", "this function will be removed. You can use 'set.cats' instead")
		set.cats(x, ...)
	}
)

