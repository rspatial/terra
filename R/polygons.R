
	
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



setMethod ('show' , 'SpatVector', 
	function(object) {
		
		cat('class       :' , class(object), '\n')
		cat('sub-class   :' , subClass(object), '\n')
		cat('geometries  : ', length(object), ' \n', sep="" ) 
		e <- as.vector(ext(object))
		cat('extent      : ' , e[1], ', ', e[2], ', ', e[3], ', ', e[4], '  (xmin, xmax, ymin, ymax)\n', sep="")
		cat('coord. ref. :' , crs(object), '\n')
		
	}
)


