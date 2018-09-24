# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# Licence GPL v3



setMethod('readStart', signature(x='SpatRaster'), 
	function(x, ...) {
		success <- try(x@ptr$readStart(), silent=TRUE)
		if (!isTRUE(success)) {
			stop("Cannot read from file")
		} else {
			invisible(success)
		}
	}
)


setMethod('readStop', signature(x='SpatRaster'), 
	function(x) {
		success <- try(x@ptr$readStop(), silent=TRUE)
		if (!isTRUE(success)) {
			stop("Cannot close an input file")
		}
		invisible(success)
	}
)


setMethod('writeStart', signature(x='SpatRaster', filename='character'), 
	function(x, filename, overwrite=FALSE, ...) {
		success <- try(x@ptr$writeStart(filename, overwrite), silent=TRUE)
		if (inherits(success, "try-error")) {
			stop("Cannot open file for writing")
		} 
		b <- x@ptr$getBlockSize()
		b$row <- b$row + 1
		b		
	}
)


setMethod('writeStop', signature(x='SpatRaster'), 
	function(x) {
		success <- try(x@ptr$writeStop(), silent=TRUE)
		if (!isTRUE(success)) {
			stop("Cannot close the output file")
		}
		invisible(success)
	} 
)


setMethod('writeValues', signature(x='SpatRaster', v='vector'), 
	function(x, v, start) {
		success <- try(x@ptr$writeValues(v, start-1), silent=TRUE)
		if (!isTRUE(success)) {
			stop("Cannot write to the output file")
		}
		invisible(success)
	}
)



setMethod('writeRaster', signature(x='SpatRaster', filename='character'), 
function(x, filename, ...) {
	if (!.hasValues(x)) {
		warning('all cell values are NA')
	}
	success <- try(x@ptr$writeRaster(filename), silent=TRUE)
	if (!isTRUE(success)) {
		stop("Cannot write to the output file")
	}
	invisible(success)
}	
)

