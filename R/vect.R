

if (!isGeneric("vect") ) { setGeneric("vect", function(x, ...) standardGeneric("vect")) }

setMethod("vect", signature(x='missing'), 
	function(...) {
		p <- methods::new('SpatLayer')
		p@ptr <- SpatLayer$new()
		.messages(p)
		return(p)
	}
)

setMethod("vect", signature(x='character'), 
	function(x, attr=NULL, crs=NA, ...) {
		p <- methods::new('SpatLayer')
		p@ptr <- SpatLayer$new()
		p@ptr$read(x)
		.messages(p)
		p
	}
)

