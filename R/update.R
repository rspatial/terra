
setMethod("update.names", signature(x="SpatRaster"),
	function(x) {
		opt <- terra:::spatOptions()
		ok <- x@ptr$update_names(opt)
		if (!ok) {
			messages(x, "update.names")
		}
		invisible(x)
	}
)

