

setMethod("set.window", signature(x="SpatRaster"),
	function(x, value)  {
		if (inherits(value, "SpatExtent")) {
			value <- value * ext(x)
			if (!(x@pntr$setWindow(value@pntr))) {
				error("window<-,SpatRaster", "could not set window")
			}
		} else if (is.null(value) || is.na(value)) {
			x@pntr$removeWindow()
		} else {
			error("window<-", "'value' should be a SpatExtent, NULL or NA")
		}
		x
	}
)


setMethod("window<-", signature(x="SpatRaster"),
	function(x, value)  {
		set.window(deepcopy(x), value)
	}
)


setMethod("window", signature(x="SpatRaster"),
	function(x)  {
		x@pntr$hasWindow()
	}
)


