

setMethod("window<-", signature(x="SpatRaster"), 
	function(x, ..., value)  {
		if (inherits(value, "SpatExtent")) {
			if (!(x@ptr$setWindow(value@ptr))) {
				stop("could not set window")
			}
		} else if (is.null(value)) {
			x@ptr$removeWindow()
		} else {
			stop("'value' should be a SpatExtent or NULL")
		}
		x
	}
)



