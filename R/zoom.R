# R function for the raster package
# Author: Robert J. Hijmans
# Date : September 2009
# Version 0.9
# Licence GPL v3


setMethod("zoom", signature(x="SpatRaster"), 
	function(x, ext=draw(), maxcell=10000, layer=1, new=TRUE, ...) {
		if (is.function(ext)) {
			ext <- ext  # force to start with drawing before creating a new graphics device
		} else {
			ext <- extent(ext)
		}
		if (new) { 
			grDevices::dev.new() 
		}
		if (nlayers(x) > 1) { 
			x <- raster(x, layer) 
		}
		plot(x, layer, col=col, maxcell=maxcell, ...) 
		return(invisible(ext))
	}
)

