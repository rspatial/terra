
setMethod("update", signature(object="SpatRaster"),
	function(object, crs=FALSE, extent=FALSE) {
		opt <- spatOptions()
		names <- FALSE
		ok <- object@pntr$update_meta(names, crs, extent, opt)
		messages(object, "update")
		invisible(object)
	}
)

