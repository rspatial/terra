# Author: Robert J. Hijmans
# Date :  April 2009
# Version 0.9
# License GPL v3


setMethod("image", signature(x="SpatRaster"), 
	function(x, y=1, maxcell=100000, xlab="", ylab="", ...)  {
		y <- as.integer(y[1])
		stopifnot(y>0 && y<=nlyr(x))
		x <- spatSample(x[[y]], maxcell, method="regular", as.raster=TRUE)
		X <- xFromCol(x, 1:ncol(x))
		Y <- yFromRow(x, nrow(x):1)
		value <- matrix(as.vector(x), nrow=nrow(x), byrow=TRUE)
		value <- t(value[nrow(value):1, ,drop=FALSE])
		if (is.null(list(...)$asp)) {
			asp <- ifelse(isLonLat(x, perhaps=TRUE, warn=FALSE), 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180), 1)
			graphics::image(x=X, y=Y, z=value, asp=asp, ...)			
			
		} else {
			graphics::image(x=X, y=Y, z=value, ...)			
		}		
	}
)

