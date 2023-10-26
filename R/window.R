

setMethod("window<-", signature(x="SpatRaster"),
	function(x, value)  {
		if (inherits(value, "SpatExtent")) {
			value <- value * ext(x)
			if (!(x@cpp$setWindow(value@cpp))) {
				error("window<-,SpatRaster", "could not set window")
			}
			#warn("window<-", "using a window is experimental")
		} else if (is.null(value) || is.na(value)) {
			x@cpp$removeWindow()
		} else {
			error("window<-", "'value' should be a SpatExtent, NULL or NA")
		}
		x
	}
)


setMethod("window", signature(x="SpatRaster"),
	function(x)  {
		x@cpp$hasWindow()
	}
)


