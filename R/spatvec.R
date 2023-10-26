


setMethod("geomtype", signature(x="SpatVector"),
	function(x){
		x@cpp$type()
	}
)
setMethod("geomtype", signature(x="SpatVectorProxy"),
	function(x){
		x@cpp$v$type()
	}
)


setMethod("datatype", signature(x="SpatVector"),
	function(x){
		x@cpp$df$get_datatypes()
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
		type <- sub("spatial", "", as.vector(tolower(class(x)[1])))
		type <- sub("dataframe", "", type)
		if (type %in% c("grid", "pixels")) type <- "raster"
		type
	}
)

setMethod("geom", signature(x="SpatVector"),
	function(x, wkt=FALSE, hex=FALSE, df=FALSE, list=FALSE, xnm="x", ynm="y"){
		if (hex) {
			x@cpp$hex()
		} else if (wkt) {
			x@cpp$getGeometryWKT()
			# or via geos with
			# x@cpp$wkt()
		} else if (list) {
			x@cpp$get_geometryList(xnm, ynm)
		} else {
			g <- x@cpp$get_geometry()
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

setMethod("crds", signature(x="SpatVector"),
	function(x, df=FALSE, list=FALSE){
		if (list) {
			gt <- geomtype(x) 
			if (gt == "lines") {
				x@cpp$linesNA()
			} else if (gt == "polygons") {
				x@cpp$polygonsList()			
			} else {
				x@cpp$coordinates()
			}
		} else {
			g <- x@cpp$coordinates()
			g <- do.call(cbind, g)
			colnames(g) <- c("x", "y")
			if (df) {
				data.frame(g)
			} else {
				g
			}
		}
	}
)


setMethod("crds", signature(x="SpatRaster"),
	function(x, df=FALSE, na.rm=TRUE, na.all=FALSE){
#		crds( as.points(x, values=(hasValues(x) || (na.rm)), na.rm=na.rm, na.all=na.all), df=df)
		opt <- spatOptions()
		out <- x@cpp$crds(na.rm, na.all, opt)
		messages(x)
		if (df) {
			out <- data.frame(out)
		} else {
			out <- do.call(cbind, out)
		}
		colnames(out) <- c("x", "y")
		out
	}
)


setMethod("dim", signature(x="SpatVector"),
	function(x){
		c(nrow(x), ncol(x))
	}
)

setMethod("dim", signature(x="SpatVectorProxy"),
	function(x){
		c(x@cpp$v$geom_count, x@cpp$v$ncol())
	}
)


as.data.frame.SpatVector <- function(x, row.names=NULL, optional=FALSE, geom=NULL, ...) {
	d <- .getSpatDF(x@cpp$df, ...)
	# fix empty names
	colnames(d)[1:ncol(x)] <- x@cpp$names
	if (!is.null(geom)) {
		geom <- match.arg(toupper(geom), c("WKT", "HEX", "XY"))
		if (geom == "XY") {
			if (!grepl("points", geomtype(x))) {
				error("as.data.frame", 'geom="XY" is only valid for point geometries')
			}
			if (nrow(d) > 0) {
				d <- cbind(d, crds(x))
			} else {
				d <- data.frame(crds(x), ...)
			}
		} else {
			g <- geom(x, wkt=geom=="WKT", hex=geom=="HEX")
			if (nrow(d) > 0) {
				d$geometry <- g
			} else {
				d <- data.frame(geometry=g, stringsAsFactors=FALSE, ...)
			}
		}
	}
	d
}
setMethod("as.data.frame", signature(x="SpatVector"), as.data.frame.SpatVector)


get.data.frame <- function(x) {
	v <- vect()
	v@cpp <- x@cpp$v
	d <- as.data.frame(v)
	d[0,,drop=FALSE]
}


as.list.SpatVector <- function(x, geom=NULL, ...) {
	as.list(as.data.frame(x, geom=geom))
}
setMethod("as.list", signature(x="SpatVector"), as.list.SpatVector)



setMethod ("expanse", "SpatVector",
	function(x, unit="m", transform=TRUE) {
		a <- x@cpp$area(unit, transform, double());
		x <- messages(x, "expanse");
		return(abs(a))
	}
)


setMethod("perim", signature(x="SpatVector"),
	function(x) {
		p <- x@cpp$length();
		x <- messages(x, "perim");
		return(p)
	}
)

setMethod("length", signature(x="SpatVector"),
	function(x) {
		x@cpp$size()
	}
)


setMethod("fillHoles", signature(x="SpatVector"),
	function(x, inverse=FALSE) {
		if (inverse) {
			x@cpp <- x@cpp$get_holes()
		} else {
			x@cpp <- x@cpp$remove_holes()
		}
		messages(x)
	}
)



#setMethod("eliminate", signature(x="SpatVector"),
#	function(x, y) {
#		x@cpp <- x@cpp$eliminate(y@cpp)
#		messages(x)
#	}
#)



setMethod("centroids", signature(x="SpatVector"),
	function(x, inside=FALSE) {
		if (inside) {
			x@cpp <- x@cpp$point_on_surface(TRUE)
		} else {
			x@cpp <- x@cpp$centroid(TRUE)
		}
		messages(x)
	}
)



setMethod("densify", signature(x="SpatVector"),
	function(x, interval, equalize=TRUE, flat=FALSE) {
		x@cpp <- x@cpp$densify(interval, equalize, flat)
		messages(x)
	}
)

setMethod("normalize.longitude", signature(x="SpatVector"),
	function(x) {
		if (nrow(x) == 0) return(deepcopy(x))
		nc <- ncol(x)
		fname <- "uuu-_123_uqq_-agg_-id123"
		x[[fname]] <- 1:nrow(x)
		d <- values(x)
		out <- x[,fname]
		out@cpp <- out@cpp$normalize_longitude()
		out <- messages(out)
		a <- aggregate(out, fname, count=FALSE)
		if (nc > 0) {
			a <- merge(a, d, by=fname)
		}
		a[[fname]] <- NULL
		a
	}
)



