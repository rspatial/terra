

setMethod("window<-", signature(x="SpatRaster"), 
	function(x, ..., value)  {
		if (inherits(value, "SpatExtent")) {
			if (!(x@ptr$setWindow(value@ptr))) {
				error("window<-,SpatRaster", "could not set window")
			}
		} else if (is.null(value)) {
			x@ptr$removeWindow()
		} else {
			error("window<-,SpatRaster", "'value' should be a SpatExtent or NULL")
		}
		x
	}
)



