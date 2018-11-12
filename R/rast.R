# Author: Robert J. Hijmans
# Date :  October 2017
# Version 1.0
# Licence GPL v3

if (!isGeneric("rast") ) { setGeneric("rast", function(x, ...) standardGeneric("rast")) }

setMethod('rast', signature(x='missing'), 
	function(x, nrows=180, ncols=360, nlyrs=1, xmin=-180, xmax=180, ymin=-90, ymax=90, crs, extent, resolution, ...) {

		if (missing(extent)) {	extent <- ext(xmin, xmax, ymin, ymax) }
		e <- as.vector(extent)
	
		if (missing(crs)) {
			if (e[1] > -360.01 & e[2] < 360.01 & e[3] > -90.01 & e[4] < 90.01) { 
				crs <- "+proj=longlat +datum=WGS84"
			} else {
				crs <- as.character(NA)
			}
		} else {
			crs <- as.character(crs)
		}
		

		r <- methods::new('SpatRaster')
		r@ptr <- SpatRaster$new(c(nrows, ncols, nlyrs), e, crs)
		
		if (!missing(resolution)) {
		#	res(r) <- resolution
			stop()
		}
		
		show_messages(r, "rast")		
	}
)


setMethod('rast', signature(x='SpatExtent'), 
	function(x, nrows=10, ncols=10, nlyrs=1, crs="", ...) {
		e <- as.vector(x)		
		r <- methods::new('SpatRaster')
		r@ptr <- SpatRaster$new(c(nrows, ncols, nlyrs), e, crs)
		show_messages(r, "rast")		
	}
)

setMethod('rast', signature(x='SpatLayer'), 
	function(x, nrows=10, ncols=10, nlyrs=1, ...) {
		rast(ext(x), nrows=nrows, ncols=ncols, nlyrs=nlyrs, crs=crs(x), ...)
	}
)



.fullFilename <- function(x, expand=FALSE) {
	x <- trimws(x)
	if (identical(basename(x), x)) {
		x <- file.path(getwd(), x)
	}
	if (expand) {
		x <- path.expand(x)
	}
	return(x)
}

setMethod('rast', signature(x='character'), 
	function(x, ...) {
		f <- .fullFilename(x)
		r <- methods::new('SpatRaster')
		r@ptr <- SpatRaster$new(f)
		show_messages(r, "rast")		
	}
)


setMethod('rast', signature(x='SpatRaster'), 
	function(x, nlyrs=nlyr(x), ...) {
		r <- methods::new('SpatRaster')
		dims <- dim(x)
		stopifnot(nlyrs > 0)
		dims[3] <- nlyrs
		r@ptr <- SpatRaster$new(dims, as.vector(ext(x)), crs(x))
		# also need the keep the names ?
		show_messages(r, "rast")		
	}
)



setMethod('rast', signature(x='matrix'), 
	function(x, ...) {
		r <- methods::new('SpatRaster')
		r@ptr <- SpatRaster$new(c(dim(x), 1), c(0, ncol(x), 0, nrow(x)), "")
		values(r) <- x
		show_messages(r, "rast")		
	}
)


setMethod('rast', signature(x='array'), 
	function(x, ...) {
		dims <- dim(x)
		if (length(dims) > 3) {
			stop("cannot handle an array with more than 3 dimensions")		
		}
		r <- methods::new('SpatRaster')
		r@ptr <- SpatRaster$new(dims, c(0, dims[2], 0, dims[1]), "")
		values(r) <- x
		show_messages(r, "rast")		
	}
)


