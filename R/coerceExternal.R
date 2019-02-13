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
# RAT tables

.fromRasterLayerBrick <- function(from) {
	f <- filename(from)
	if (f != "") {
		r <- rast(f)
		if (from@file@NAchanged) {
			warning("changed NA value ignored")
		}
		return(r)
	} else {
		r <- rast(	nrows=nrow(from), 
					ncols=ncol(from),
					nlyrs=nlayers(from),
					crs=crs(from),
					extent=extent(from))
		if (hasValues(from)) {
			values(r) <- values(from)			
		}	
		names(r)  <- names(from)
	}
	return(r)
}

.fromRasterStack <- function(from) {
	x <- from[[1]]
	n <- nbands(x)
	if ((n > 1) & (n == nlayers(from))) {
		ff <- lapply(1:nlayers(from), function(i) { filename(from[[i]]) })
		if (length(unique(ff)) == 1) {
			r <- rast(filename(x))
			return(r)
		}	
	} 
	s <- lapply(1:nlayers(from), function(i) {
		x <- from[[i]]
		.fromRasterLayerBrick(x)[[bandnr(x)]]
	})
	do.call(c, s)
}


setAs("Raster", "SpatRaster", 
	function(from) {
		if (inherits(from, "RasterLayer") | inherits(from, "RasterBrick")) { 
			.fromRasterLayerBrick(from)			
		} else {
			.fromRasterStack(from)
		}
	}
)

