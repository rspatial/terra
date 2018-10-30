
	
setMethod ('length' , 'SpatVector', 
	function(x) {
		x@ptr$size()
	}
)

setMethod('ext', signature(x='SpatVector'), 
	function(x, ...){ 
		e <- methods::new('SpatExtent')
		e@ptr <- x@ptr$extent
		return(e)
	}
)	


if (!isGeneric("subClass")) {setGeneric("subClass", function(x,...) standardGeneric("subClass"))}	

setMethod('subClass', signature(x='SpatVector'), 
	function(x, ...){ 
		a <- as.vector(class(x@ptr))
		tolower(gsub("Rcpp_Spat", "", a))
	}
)	


