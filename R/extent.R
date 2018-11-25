# Author: Robert J. Hijmans
# Date :  October 2017
# Version 0.9
# License GPL v3
	

setMethod('ext', signature(x='missing'), 
	function(x, ...){ 
		e <- methods::new('SpatExtent')
		e@ptr <- SpatExtent$new()
		return(e)
	}
)	
	
setMethod('ext', signature(x='numeric'), 
	function(x, ...){ 
		dots <- unlist(list(...))
		x <- c(x, dots)
		if (length(x) < 4) {
			stop('insufficient number of elements (should be 4)')
		}
		if (length(x) > 4) {
			warning('more elements than expected (should be 4)')
		}
		names(x) <- NULL
		e <- methods::new('SpatExtent')
		e@ptr <- SpatExtent$new(x[1], x[2], x[3], x[4])
		if (validObject(e)) return(e)
	}	
)


setMethod('ext', signature(x='SpatRaster'), 
	function(x, ...){ 
		e <- methods::new('SpatExtent')
		e@ptr <- x@ptr$extent
		return(e)
	}
)	

setMethod("ext<-", signature('SpatRaster', 'SpatExtent'), 
	function(x, value) {
	    stop("not yet implemented")
		#y@ptr <- y@ptr$extent(value)
		#show_messages(y)
	}
)



setMethod('ext', signature(x='SpatVector'), 
	function(x, ...){ 
		e <- methods::new('SpatExtent')
		e@ptr <- x@ptr$extent()
		return(e)
	}
)	
