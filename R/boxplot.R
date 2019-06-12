# Author: Robert J. Hijmans 
# Date :  November 2010
# Version 1.0
# Licence GPL v3
 

if (!isGeneric("boxplot")) {
	setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))
}

setMethod("boxplot", signature(x="SpatRaster"), 
	function(x, maxpixels=100000, ...) {
		cn <- names(x)
		if ( ncell(x) > maxpixels) {
			warning("taking a sample of ", maxpixels, " cells")
			x <- sampleRegular(x, maxpixels)
		} 
		x <- values(x)
		colnames(x) <- cn
		boxplot(x, ...)
	}
)


if (!isGeneric("barplot")) {
	setGeneric("barplot", function(height,...) standardGeneric("barplot"))
}	

setMethod("barplot", "SpatRaster", 
	function(height, maxpixels=1000000, digits=0, breaks=NULL, col=grDevices::rainbow, ...) {
		
		x <- values(sampleRegular(height[[1]], maxpixels))
		adj <- length(x) / ncell(height)
		if (adj < 1) {
			warning('a sample of ', round(100*adj, 1), '% of the raster cells were used to estimate frequencies')
		}

		if (!is.null(digits)) {
			x <- round(x, digits)
		}
		
		if (!is.null(breaks)) {
			x <- cut(x, breaks)
		}
		
		x <- table(x) / adj
		if (is.function(col)) {
			col <- col(length(x))
		}
		barplot(x, col=col, ...)
	}
)
