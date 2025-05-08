

setMethod("set.window", signature(x="SpatRaster"),
	function(x, value)  {
		try(value <- ext(value), silent=TRUE)
		if (inherits(value, "SpatExtent")) {
			window(x) <- NULL
			value <- intersect(ext(x), value)
			if (is.null(value)) {
				error("window<-,SpatRaster", "cannot set a window that does not overlap with x")			
			}
			if (!(x@pntr$setWindow(value@pntr))) {
				error("window<-,SpatRaster", "cannot set this window")
			}
		} else if (is.null(value) || isTRUE(is.na(value))) {
			x@pntr$removeWindow()
		} else {
			error("window<-", "'value' should be or have a SpatExtent, or be NULL or NA")
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


