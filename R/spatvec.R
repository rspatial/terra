
if (!isGeneric("geomtype")) {setGeneric("geomtype", function(x,...) standardGeneric("geomtype"))}	
if (!isGeneric("geometry")) {setGeneric("geometry", function(x,...) standardGeneric("geometry"))}	


setMethod ('length' , 'SpatLayer', 
	function(x) {
		x@ptr$size()
	}
)

setMethod('geomtype', signature(x='SpatLayer'), 
	function(x, ...){ 
		x@ptr$type()
	}
)	

setMethod('geometry', signature(x='SpatLayer'), 
	function(x, ...){ 
		x@ptr$getGeometry()
	}
)	

setMethod('as.data.frame', signature(x='SpatLayer'), 
	function(x, ...) {
		d <- data.frame(x@ptr$getAttributes(), stringsAsFactors=FALSE)
		names(d) <- x@ptr$names()
		d
	}
)
	
	