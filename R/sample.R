

setMethod("sampleRegular", signature(x="SpatRaster", size="numeric"), 
	function(x, size, type="regular", ...) { 
		size <- max(1, min(size(x), size))
		x@ptr <- x@ptr$sampleRegular(size)
		show_messages(x, "sampleRegular")		
	}
)


