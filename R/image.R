# Author: Robert J. Hijmans
# Date :  April 2019
# Version 1.0
# License GPL v3

.plot_image <- function(x, y=1, xlab="", ylab="", asp=NULL, ...) {
	X <- xFromCol(x, 1:ncol(x))
	Y <- yFromRow(x, nrow(x):1)
	value <- matrix(as.vector(x), nrow=nrow(x), byrow=TRUE)
	value <- t(value[nrow(value):1, ,drop=FALSE])
	if (is.null(asp)) {
		asp <- ifelse(is.lonlat(x, perhaps=TRUE, warn=FALSE), 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180), 1)
		graphics::image(x=X, y=Y, z=value, asp=asp, xlab=xlab, ylab=ylab, ...)
	} else {
		graphics::image(x=X, y=Y, z=value, xlab=xlab, ylab=ylab, ...)
	}
}

setMethod("image", signature(x="SpatRaster"),
	function(x, y=1, maxcell=500000, ...)  {
		y <- as.integer(y[1])
		stopifnot(y > 0 && y <= nlyr(x))
		x <- spatSample(x[[y]], maxcell, method="regular", as.raster=TRUE, warn=FALSE)
		.plot_image(x, ...)
	}
)

