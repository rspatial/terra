# Author: Robert J. Hijmans
# Date : November 2009
# Version 1.0
# Licence GPL v3


setMethod('rasterize', signature(x='SpatPolygons', y='SpatRaster'), 
function(x, y, background=NA, filename="", ...) { 
	overwrite <- .overwrite(...)
	r <- methods::new('SpatRaster')
	r@ptr <- y@ptr$rasterizePolygons(x@ptr, background, filename[1], overwrite)
	return(r)
}
)
