
setMethod("update", signature(object="SpatRaster"),
	function(object, names=FALSE) {
		if (names) {
			opt <- terra:::spatOptions()
			ok <- object@ptr$update_names(opt)
			if (!ok) {
				messages(object, "update")
			}
		}
		invisible(object)
	}
)

