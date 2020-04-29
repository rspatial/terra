# R function for the raster package
# Author: Robert J. Hijmans
# Date : September 2009
# Version 0.9
# Licence GPL v3


setMethod("zoom", signature(x="SpatRaster"), 
	function(x, e=draw(), maxcell=10000, layer=1, new=TRUE, ...) {
		if (is.function(ext)) {
			ext <- e  # force to start with drawing before creating a new graphics device
		} else {
			if (!inherits(e, "SpatExtent")) {
				e <- ext(e)
			}
		}
		if (new) { 
			grDevices::dev.new() 
		}
		x <- crop(x[[layer]], e)
		plot(x, maxcell=maxcell, ...) 
		return(invisible(e))
	}
)

