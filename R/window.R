

setMethod("window<-", signature(x="SpatRaster"), 
	function(x, ..., value)  {
		if (inherits(value, "SpatExtent")) {
			value <- value * ext(x)
			if (!(x@ptr$setWindow(value@ptr))) {
				error("window<-,SpatRaster", "could not set window")
			}
			warn("window<-", "using a window is experimental. User beware")
		} else if (is.null(value	)) {
			x@ptr$removeWindow()
		} else {
			error("window<-,SpatRaster", "'value' should be a SpatExtent or NULL")
		}
		x
	}
)


setMethod("window", signature(x="SpatRaster"), 
	function(x, ...)  {
		x@ptr$hasWindow()
	}
)


