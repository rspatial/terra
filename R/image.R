# Author: Robert J. Hijmans
# Date :  April 2009
# Version 0.9
# Licence GPL v3

if (!isGeneric("image")) {
	setGeneric("image", function(x,...)
		standardGeneric("image"))
}	

setMethod("image", signature(x='SpatRaster'), 
	function(x, maxpixels=100000, xlab='', ylab='', useRaster=TRUE, ...)  {
		x <- sampleRegular(x[[1]], maxpixels)
		X <- xFromCol(x, 1:ncol(x))
		Y <- yFromRow(x, nrow(x):1)
		value <- matrix(as.vector(x), nrow=nrow(x), byrow=TRUE)
		value <- t(value[nrow(value):1, ,drop=FALSE])
		image(x=X, y=Y, z=value, useRaster=useRaster, ...)			
	}
)

