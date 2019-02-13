# Author: Robert J. Hijmans 
# Date : February 2019
# Version 1.0
# License GPL v3

# todo:
# for ncdf files (not yet natively supported in terra)
# check the variable to be used
# 
# check z values, other attributes such as NAvalue that may have been
# changed after creation of object from file

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
		names(r)  <- names(from)
	}
	return(r)
}

setAs("Raster", "SpatRaster", 
	function(from) {
		if (inherits(from, "RasterLayer") | inherits(from, "RasterBrick")) { 
			.fromRasterLayerBrick(from)			
		} else {
			if (raster::canProcessInMemory(from)) {
				# here we could first check if all bands are used
				# in the right order, to create a more efficient object
				# but for now:
				s <- lapply(1:nlayers(from), function(i) {
					x <- from[[i]]
					.fromRasterLayerBrick(x)[[bandnr(x)]]
				})
				do.call(c, s)
			}
		}
	}
)

