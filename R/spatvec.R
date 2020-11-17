

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


setMethod("is.lines", signature(x="SpatVector"), 
	function(x, ...) {
		geomtype(x) == "lines"
	}
)

setMethod("is.polygons", signature(x="SpatVector"), 
	function(x, ...) {
		geomtype(x) == "polygons"
	}
)
setMethod("is.points", signature(x="SpatVector"), 
	function(x, ...) {
		geomtype(x) == "points"
	}
)


setMethod("geomtype", signature(x="Spatial"), 
	function(x, ...){ 
		type <- sub("spatial", "", as.vector(tolower(class(x))))
		type <- sub("dataframe", "", type)
		if (type %in% c("grid", "pixels")) type <- "raster"
		type
	}
)	

setMethod("geom", signature(x="SpatVector"), 
	function(x, wkt=FALSE, ...){
		if (wkt) {
			x@ptr$getGeometryWKT()
		} else {
			x@ptr$getGeometry()
		}
	}
)	


setMethod("dim", signature(x="SpatVector"), 
	function(x){ 
		c(nrow(x), ncol(x))
	}
)	

setMethod("as.data.frame", signature(x="SpatVector"), 
	function(x, geom=FALSE, ...) {
		d <- data.frame(x@ptr$getDF(), check.names=FALSE, fix.empty.names=FALSE, stringsAsFactors=FALSE)
		colnames(d) <- x@ptr$names
		if (geom) {
			d$geometry <- geom(x, TRUE)
		}
		d
	}
)

setMethod("as.list", signature(x="SpatVector"), 
	function(x, geom=FALSE, ...) {
		as.list(as.data.frame(x, geom=geom, ...))
	}
)
	
	

setMethod("area", signature(x="SpatVector"), 
	function(x, ...) {
		a <- x@ptr$area();
		x <- show_messages(x, "area");
		return(a)
	}
)	

setMethod("perimeter", signature(x="SpatVector"), 
	function(x) {
		a <- x@ptr$length();
		x <- show_messages(x, "length");
		return(a)
	}
)	

setMethod("length", signature(x="SpatVector"), 
	function(x) {
		size(x)
	}
)	



setMethod("fill", signature(x="SpatVector"), 
	function(x, ...) {
		x@ptr <- x@ptr$remove_holes()
		x
	}
)	
