
setMethod("update.names", signature(x="SpatRaster"),
	function(x) {
		opt <- terra:::spatOptions()
		ok <- x@ptr$update_names(opt)
		if (!ok) {
			error("update.names", "update failed")
		}
		invisible(x)
	}
)
