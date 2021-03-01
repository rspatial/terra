

setMethod ("size" , "SpatVector", 
	function(x) {
		x@ptr$size()
	}
)

setMethod("geomtype", signature(x="SpatVector"), 
	function(x){ 
		x@ptr$type()
	}
)

setMethod("datatype", signature(x="SpatVector"), 
	function(x){ 
		x@ptr$df$get_datatypes()
	}
)


setMethod("is.lines", signature(x="SpatVector"), 
	function(x) {
		geomtype(x) == "lines"
	}
)

setMethod("is.polygons", signature(x="SpatVector"), 
	function(x) {
		geomtype(x) == "polygons"
	}
)
setMethod("is.points", signature(x="SpatVector"), 
	function(x) {
		grepl("points", geomtype(x))
	}
)


setMethod("geomtype", signature(x="Spatial"), 
	function(x){ 
		type <- sub("spatial", "", as.vector(tolower(class(x))))
		type <- sub("dataframe", "", type)
		if (type %in% c("grid", "pixels")) type <- "raster"
		type
	}
)

setMethod("geom", signature(x="SpatVector"), 
	function(x, wkt=FALSE, df=FALSE){
		if (wkt) {
			x@ptr$getGeometryWKT()
		} else {
			g <- x@ptr$get_geometry()
			g <- do.call(cbind, g)
			colnames(g) <- c("geom", "part", "x", "y", "hole")[1:ncol(g)]
			if (df) {
				data.frame(g)
			} else {
				g
			}
		}
	}
)

setMethod("coords", signature(x="SpatVector"), 
	function(x, df=FALSE){
		g <- x@ptr$coordinates()
		g <- do.call(cbind, g)
		colnames(g) <- c("x", "y")
		if (df) {
			data.frame(g)
		} else {
			g
		}
	}
)

setMethod("coords", signature(x="SpatRaster"), 
	function(x, df=FALSE){
		x <- as.points(x)
		coordinates(x)
	}
)


setMethod("dim", signature(x="SpatVector"), 
	function(x){ 
		c(nrow(x), ncol(x))
	}
)

setMethod("as.data.frame", signature(x="SpatVector"), 
	function(x, geom=FALSE) {
		d <- data.frame(x@ptr$getDF(), check.names=FALSE, fix.empty.names=FALSE, stringsAsFactors=FALSE)
		colnames(d) <- x@ptr$names
		if (geom) {
			g <- geom(x, wkt=TRUE)
			if (nrow(d) > 0) {
				d$geometry <- g
			} else {
				d <- data.frame(geometry=g, stringsAsFactors=FALSE)
			}
		}
		d
	}
)

setMethod("as.list", signature(x="SpatVector"), 
	function(x, geom=FALSE) {
		as.list(as.data.frame(x, geom=geom))
	}
)



setMethod("area", signature(x="SpatVector"), 
	function(x) {
		a <- x@ptr$area();
		x <- messages(x, "area");
		return(a)
	}
)

setMethod("perimeter", signature(x="SpatVector"), 
	function(x) {
		p <- x@ptr$length();
		x <- messages(x, "length");
		return(p)
	}
)

setMethod("length", signature(x="SpatVector"), 
	function(x) {
		size(x)
	}
)



setMethod("fill", signature(x="SpatVector"), 
	function(x, inverse=FALSE) {
		if (inverse) {
			x@ptr <- x@ptr$get_holes()
		} else {
			x@ptr <- x@ptr$remove_holes()
		}
		messages(x)
	}
)



setMethod("centroids", signature(x="SpatVector"), 
	function(x) {
		x@ptr <- x@ptr$centroid()
		messages(x)
	}
)

