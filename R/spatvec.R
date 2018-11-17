

setMethod ("size" , "SpatVector", 
	function(x, ...) {
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
	

setMethod("area", signature(x="SpatVector"), 
	function(x, ...) {
		a <- x@ptr$area();
		x <- show_messages(x, "area");
		return(a)
	}
)	

setMethod("length", signature(x="SpatVector"), 
	function(x) {
		a <- x@ptr$length();
		x <- show_messages(x, "length");
		return(a)
	}
)	
