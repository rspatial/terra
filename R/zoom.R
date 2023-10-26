# Author: Robert J. Hijmans
# Date : September 2009
# Version 0.9
# License GPL v3


setMethod("zoom", signature(x="SpatRaster"),
	function(x, e=draw(), maxcell=100000, layer=1, new=FALSE, ...) {
		if (grDevices::dev.cur() == 1) {
			if (!is.null(RGB(x))) {
				plot(x, maxcell=maxcell, ...)
			} else {
				plot(x, layer, maxcell=maxcell, ...)
			}
		}
		if (is.function(e)) {
		# force to start with drawing before creating a new graphics device
			e <- e
		} else if (!inherits(e, "SpatExtent")) {
			e <- ext(e)
		}
		if (new) {
			grDevices::dev.new()
		}
		if (!is.null(RGB(x))) {
			x <- crop(x, e)
		} else {
			x <- crop(x[[layer]], e)
		}
		plot(x, maxcell=maxcell, ...)
		return(invisible(e))
	}
)


setMethod("zoom", signature(x="SpatVector"),
	function(x, e=draw(), new=FALSE, ...) {

		if (grDevices::dev.cur() == 1) {
			plot(x, ...)
		}
		if (is.function(e)) {
			e <- e
		} else if (!inherits(e, "SpatExtent")) {
			e <- ext(e)
		}
		if (new) {
			grDevices::dev.new()
		}
		x <- crop(x, e)
		plot(x, ...)
		return(invisible(e))
	}
)
