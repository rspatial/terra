# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3



setMethod('readStart', signature(x='SpatRaster'), 
	function(x, ...) {
		success <- x@ptr$readStart()
		show_messages(x, "readStart")		
		invisible(success)
	}
)


setMethod('readStop', signature(x='SpatRaster'), 
	function(x) {
		success <- x@ptr$readStop()
		show_messages(x, "readStop")		
		invisible(success)
	}
)


setMethod('writeStart', signature(x='SpatRaster', filename='character'), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename[1], overwrite[1], wopt)
		x@ptr$writeStart(opt)
		show_messages(x, "writeStart")		
		b <- x@ptr$getBlockSize(4)
		b$row <- b$row + 1
		b		
	}
)


setMethod('writeStop', signature(x='SpatRaster'), 
	function(x) {
		success <- x@ptr$writeStop()
		show_messages(x, "writeStop")		
		invisible(success)
	} 
)


setMethod('writeValues', signature(x='SpatRaster', v='vector'), 
	function(x, v, start) {
		success <- x@ptr$writeValues(v, start-1)
		show_messages(x, "writeValues")
		invisible(success)
	}
)



setMethod('writeRaster', signature(x='SpatRaster', filename='character'), 
function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
	opt <- .runOptions(filename[1], overwrite[1],wopt)
	success <- x@ptr$writeRaster(opt)
	show_messages(x, "writeRaster")
	invisible(rast(filename))
}
)

