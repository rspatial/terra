
setMethod("update", signature(object="SpatRaster"),
	function(object, names=FALSE, crs=FALSE, extent=FALSE) {
		opt <- spatOptions()
		ok <- object@pnt$update_meta(names, crs, extent, opt)
		messages(object, "update")
		invisible(object)
	}
)

