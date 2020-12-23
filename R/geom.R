

roundtrip <- function(x, coll=FALSE) {
	if (coll) {
		p <- methods::new("SpatVectorCollection")
		p@ptr <- x@ptr$bienvenue()	
		return(p)
	} else {
		x@ptr <- x@ptr$allerretour()
		return(x)
	}
}

setMethod("is.valid", signature(x="SpatVector"), 
	function(x, messages=FALSE, ...) {
		if (messages) {
			x@ptr$geos_isvalid_msg()
		} else {
			x@ptr$geos_isvalid()
		}
	}
)



setMethod("intersect", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x@ptr <- x@ptr$intersect(y@ptr)
		messages(x)
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatExtent"), 
	function(x, y) {
		x * y
	}
)

setMethod("intersect", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y) {
		y <- as.polygons(y)
		intersect(x, y)
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatVector"), 
	function(x, y) {
		y <- ext(y)
		x * y
	}
)


setMethod("intersects", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		out <- x@ptr$intersects(y@ptr)
		x <- messages(x)
		matrix(out, nrow=nrow(x), byrow=TRUE)
	}
)


setMethod("buffer", signature(x="SpatVector"), 
	function(x, width, quadsegs=10, capstyle="round", ...) {
		if (geomtype(x) == "points") {
			x@ptr <- x@ptr$buffer(width, quadsegs, 1)
		} else {
			if (is.lonlat(x)) {
				warn("buffer", "lon/lat data are treated as planar")
			}
			quadsegs <- max(3, quadsegs)
			x@ptr <- x@ptr$buffer2(width, quadsegs, 1)		
		}
		messages(x, "buffer")
	}
)


setMethod("crop", signature(x="SpatVector", y="ANY"), 
	function(x, y, ...) {
		if (!inherits(y, "SpatExtent")) {
			y <- try(ext(y), silent=TRUE)
			if (inherits(y, "try-error")) {
				stop("y does not have a SpatExtent")
			}
		}
		x@ptr <- x@ptr$crop(y@ptr)
		messages(x)
	}
)


setMethod("disaggregate", signature(x="SpatVector"), 
	function(x, ...) {
		x@ptr <- x@ptr$disaggregate()
		messages(x, "disaggregate")
	}
)



setMethod("voronoi", signature(x="SpatVector"), 
	function(x, bnd=NULL, tolerance=0, as.lines=FALSE, ...) {
		if (is.null(bnd)) {
			bnd <- vect()
		} 
		if (inherits(bnd, "SpatExtent")) {
			bnd <- as.polygons(bnd)
		}
		x@ptr <- x@ptr$voronoi(bnd@ptr, tolerance, as.lines)
		messages(x, "voronoi")
	}
)

