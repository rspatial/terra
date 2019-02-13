# Author: Robert J. Hijmans 
# Date : October 2018
# Version 1.0
# License GPL v3


.fromRasterLayerBrick <- function(from) {
	f <- filename(from)
	if (f != "") {
		r <- rast(f)
		return(r)
	} else {
		r <- rast(	nrows=nrow(from), 
					ncols=ncol(from),
					nlyrs=nlayers(from),
					crs=crs(from),
					extent=extent(from))
		values(r) <- values(from)			
		names(r) <- names(from)
	}
	return(r)
}

setAs("Raster", "SpatRaster", 
	function(from) {
		if (inherits(from, "RasterLayer") | inherits(from, "RasterBrick")) { 
			.fromRasterLayerBrick(from)			
		} else {
			if (raster::canProcessInMemory(from)) {
				s <- lapply(1:nlayers(from), function(i) .fromRasterLayerBrick(from[[i]]))
				do.call(c, s)
			}
		}
	}
)

