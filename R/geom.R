

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
			r <- x@ptr$geos_isvalid_msg()
			d <- data.frame(matrix(r, ncol=2, byrow=TRUE))
			d[,1] = d[,1] == "\001"
			colnames(d) <- c("valid", "reason")
			d
		} else {
			x@ptr$geos_isvalid()
		}
	}
)

setMethod("cover", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, identity=FALSE, ...) {
		x@ptr <- x@ptr$cover(y@ptr, identity[1])
		messages(x, "cover")
	}
)


setMethod("symdif", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, ...) {
		x@ptr <- x@ptr$symdif(y@ptr)
		messages(x, "symdif")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, ...) {
		x@ptr <- x@ptr$erase(y@ptr)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y, ...) {
		y <- as.polygons(y)
		x@ptr <- x@ptr$erase(y@ptr)
		messages(x, "erase")
	}
)


setMethod("union", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x@ptr <- x@ptr$union(y@ptr)
		messages(x, "union")
	}
)

setMethod("union", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y) {
		y <- as.vector(y)
		x@ptr <- x@ptr$union(y@ptr)
		messages(x, "union")
	}
)

setMethod("union", signature(x="SpatExtent", y="SpatExtent"), 
	function(x, y) {
		x + y
	}
)


setMethod("intersect", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x@ptr <- x@ptr$intersect(y@ptr)
		messages(x, "intersect")
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
		x@ptr <- x@ptr$crop_ext(y@ptr)
		messages(x, "crop")
	}
)

setMethod("crop", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, ...) {
		x@ptr <- x@ptr$crop_vct(y@ptr)
		messages(x, "crop")
	}
)

setMethod("convexhull", signature(x="SpatVector"), 
	function(x, ...) {
		x@ptr <- x@ptr$chull()
		messages(x, "convexhull")
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
		} else if (inherits(bnd, "SpatExtent")) {
			bnd <- as.polygons(bnd)
		}
		x@ptr <- x@ptr$voronoi(bnd@ptr, tolerance, as.lines)
		messages(x, "voronoi")
	}
)


setMethod("delauny", signature(x="SpatVector"), 
	function(x, tolerance=0, as.lines=FALSE, ...) {
		x@ptr <- x@ptr$delauny(tolerance, as.lines)
		messages(x, "delauny")
	}
)

