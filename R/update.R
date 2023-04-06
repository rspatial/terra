
setMethod("update", signature(object="SpatRaster"),
	function(object, names=FALSE, crs=FALSE, extent=FALSE) {
		opt <- spatOptions()
		ok <- object@ptr$update_meta(names, crs, extent, opt)
		messages(object, "update")
		invisible(object)
	}
)

