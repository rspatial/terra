# Author: Robert J. Hijmans
# Date :  October 2017
# Version 0.9
# License GPL v3
	

setMethod("ext", signature(x="missing"), 
	function(x, ...){ 
		e <- methods::new("SpatExtent")
		e@ptr <- SpatExtent$new()
		return(e)
	}
)	
	
setMethod("ext", signature(x="numeric"), 
	function(x, ...){ 
		dots <- unlist(list(...))
		x <- c(x, dots)
		if (length(x) < 4) {
			stop("insufficient number of elements (should be 4)")
		}
		if (length(x) > 4) {
			warning("more elements than expected (should be 4)")
		}
		names(x) <- NULL
		e <- methods::new("SpatExtent")
		e@ptr <- SpatExtent$new(x[1], x[2], x[3], x[4])
		if (methods::validObject(e)) return(e)
	}	
)


setMethod("ext", signature(x="SpatRaster"), 
	function(x, ...){ 
		e <- methods::new("SpatExtent")
		e@ptr <- x@ptr$extent
		return(e)
	}
)	


setMethod("ext", signature(x="SpatDataSet"), 
	function(x, ...){ 
		e <- methods::new("SpatExtent")
		e@ptr <- x[1]@ptr$extent
		return(e)
	}
)	


setMethod("ext<-", signature("SpatRaster", "SpatExtent"), 
	function(x, value) {
		x@ptr$extent <- value@ptr
		show_messages(x)
	}
)


setMethod("ext<-", signature("SpatRaster", "numeric"), 
	function(x, value) {
		stopifnot(length(value) == 4)
		e <- ext(value[1], value[2], value[3], value[4])
		if (!e@ptr$valid) {
			stop("not a valid extent specification")
		}
		x@ptr$extent <- e@ptr
		show_messages(x)
	}
)



setMethod("ext", signature(x="SpatVector"), 
	function(x, ...) { 
		e <- methods::new("SpatExtent")
		e@ptr <- x@ptr$extent()
		return(e)
	}
)	


setMethod("ext", signature(x="Extent"), 
	function(x, ...) {
		ext(as.vector(x))
	}
)	

setMethod("ext", signature(x="Raster"), 
	function(x, ...) { 	
		ext(x@extent)
	}
)	

setMethod("ext", signature(x="Spatial"), 
	function(x, ...) { 	
		ext(as.vector(t(x@bbox)))
	}
)	



setMethod("xmin", signature(x="SpatExtent"), 
	function(x){ 
		x@ptr$vector[1]
	}
)	
setMethod("xmax", signature(x="SpatExtent"), 
	function(x){ 
		x@ptr$vector[2]
	}
)	
setMethod("ymin", signature(x="SpatExtent"), 
	function(x){ 
		x@ptr$vector[3]
	}
)	
setMethod("ymax", signature(x="SpatExtent"), 
	function(x){ 
		x@ptr$vector[4]
	}
)	


setMethod("xmin<-", signature("SpatExtent", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(x)
		v[1] <- value
		ext(v)
	}
)	
setMethod("xmax<-", signature("SpatExtent", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(x)
		v[2] <- value
		ext(v)
	}
)	
setMethod("ymin<-", signature("SpatExtent", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(x)
		v[3] <- value
		ext(v)
	}
)	
setMethod("ymax<-", signature("SpatExtent", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(x)
		v[4] <- value
		ext(v)
	}
)	


setMethod("xmin", signature(x="SpatRaster"), 
	function(x){ 
		xmin(ext(x))
	}
)	
setMethod("xmax", signature(x="SpatRaster"), 
	function(x){ 
		xmax(ext(x))
	}
)	
setMethod("ymin", signature(x="SpatRaster"), 
	function(x){ 
		ymin(ext(x))
	}
)	
setMethod("ymax", signature(x="SpatRaster"), 
	function(x){ 
		ymax(ext(x))
	}
)	


setMethod("xmin<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(ext(x))
		v[1] <- value
		ext(x) <- ext(v)
		x
	}
)	


setMethod("xmax<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(ext(x))
		v[2] <- value
		ext(x) <- ext(v)
		x
	}
)	
setMethod("ymin<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(ext(x))
		v[3] <- value
		ext(x) <- ext(v)
		x
	}
)	

setMethod("ymax<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value){ 
		v <- as.vector(ext(x))
		v[4] <- value
		ext(x) <- ext(v)
		x
	}
)	



setMethod("xmin", signature(x="SpatVector"), 
	function(x){ 
		xmin(ext(x))
	}
)	
setMethod("xmax", signature(x="SpatVector"), 
	function(x){ 
		xmax(ext(x))
	}
)	
setMethod("ymin", signature(x="SpatVector"), 
	function(x){ 
		ymin(ext(x))
	}
)	
setMethod("ymax", signature(x="SpatVector"), 
	function(x){ 
		ymax(ext(x))
	}
)	

.ext2bb <- function(e) {
	matrix(as.vector(e), ncol=2, byrow=TRUE)
}

setMethod("bbox", signature(obj="SpatRaster"), 
	function(obj){ 
		.ext2bb(ext(obj))
	}
)	

setMethod("bbox", signature(obj="SpatVector"), 
	function(obj){ 
		.ext2bb(ext(obj))
	}
)	
