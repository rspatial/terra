# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# Licence GPL v3



setMethod('readStart', signature(x='SpatRaster'), 
	function(x, ...) {
		success <- x@ptr$readStart()
		.messages(x, "readStart")		
		invisible(success)
	}
)


setMethod('readStop', signature(x='SpatRaster'), 
	function(x) {
		success <- x@ptr$readStop()
		.messages(x, "readStop")		
		invisible(success)
	}
)


setMethod('writeStart', signature(x='SpatRaster', filename='character'), 
	function(x, filename, overwrite=FALSE, ...) {
		x@ptr$writeStart(filename, overwrite)
		.messages(x, "writeStart")		
		b <- x@ptr$getBlockSize(4)
		b$row <- b$row + 1
		b		
	}
)


setMethod('writeStop', signature(x='SpatRaster'), 
	function(x) {
		succes <- x@ptr$writeStop()
		.messages(x, "writeStop")		
		invisible(success)
	} 
)


setMethod('writeValues', signature(x='SpatRaster', v='vector'), 
	function(x, v, start) {
		success <- x@ptr$writeValues(v, start-1)
		.messages(x, "writeValues")
		invisible(success)
	}
)



setMethod('writeRaster', signature(x='SpatRaster', filename='character'), 
function(x, filename, overwrite=FALSE, ...) {
	success <- x@ptr$writeRaster(filename, overwrite=overwrite)
	.messages(x, "writeRaster")
	invisible(rast(filename))
}
)

