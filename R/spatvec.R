

setMethod ("length" , "SpatVector", 
	function(x) {
		x@ptr$size()
	}
)

setMethod("geomtype", signature(x="SpatVector"), 
	function(x, ...){ 
		x@ptr$type()
	}
)	

setMethod("geom", signature(x="SpatVector"), 
	function(x, ...){ 
		x@ptr$getGeometry()
	}
)	

setMethod("as.data.frame", signature(x="SpatVector"), 
	function(x, ...) {
		d <- data.frame(x@ptr$getAttributes(), stringsAsFactors=FALSE)
		names(d) <- x@ptr$names()
		d
	}
)
	
	