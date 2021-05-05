

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
	function(x, messages=FALSE) {
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
	function(x, y, identity=FALSE) {
		x@ptr <- x@ptr$cover(y@ptr, identity[1])
		messages(x, "cover")
	}
)


setMethod("symdif", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x@ptr <- x@ptr$symdif(y@ptr)
		messages(x, "symdif")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x@ptr <- x@ptr$erase(y@ptr)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y) {
		y <- as.polygons(y)
		x@ptr <- x@ptr$erase(y@ptr)
		messages(x, "erase")
	}
)



setMethod("union", signature(x="SpatVector", y="missing"), 
	function(x, y) {
		x@ptr <- x@ptr$union_self()
		messages(x, "union")
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
		x@ptr$intersect(y@ptr)
		x
	}
)

setMethod("intersect", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y) {
		x@ptr <- x@ptr$crop_ext(y@ptr)
		x
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatVector"), 
	function(x, y) {
		y <- ext(y)
		x * y
	}
)

#setMethod("intersect", signature(x="SpatRaster", y="SpatRaster"),
#	function(x, y) {
#		a <- crop(x, y)
#		b <- crop(y, x)
#		c(a, b)
#	}
#)

setMethod("buffer", signature(x="SpatVector"), 
	function(x, width, quadsegs=10) {
		x@ptr <- x@ptr$buffer(width, quadsegs)
		messages(x, "buffer")
	}
)


setMethod("crop", signature(x="SpatVector", y="ANY"), 
	function(x, y) {
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
	function(x, y) {
		if (size(y) > 1) {
			y <- aggregate(y)
		}
		x@ptr <- x@ptr$crop_vct(y@ptr)
		messages(x, "crop")
	}
)

setMethod("convhull", signature(x="SpatVector"), 
	function(x) {
		x@ptr <- x@ptr$chull()
		messages(x, "convhull")
	}
)


setMethod("disaggregate", signature(x="SpatVector"), 
	function(x) {
		x@ptr <- x@ptr$disaggregate()
		messages(x, "disaggregate")
	}
)


setMethod("voronoi", signature(x="SpatVector"), 
	function(x, bnd=NULL, tolerance=0, as.lines=FALSE) {
		if (is.null(bnd)) {
			bnd <- vect()
		} else {
			bnd <- as.polygons(ext(bnd))
		}
		x@ptr <- x@ptr$voronoi(bnd@ptr, tolerance, as.lines)
		messages(x, "voronoi")
	}
)


setMethod("delauny", signature(x="SpatVector"), 
	function(x, tolerance=0, as.lines=FALSE) {
		x@ptr <- x@ptr$delauny(tolerance, as.lines)
		messages(x, "delauny")
	}
)


setMethod("flip", signature(x="SpatVector"), 
	function(x, direction="vertical") {
		d <- match.arg(direction, c("vertical", "horizontal")) 
		x@ptr <- x@ptr$flip(d == "vertical")
		messages(x, "flip")
	}
)



setMethod("spin", signature(x="SpatVector"), 
	function(x, angle, x0, y0) { 
		e <- as.vector(ext(x))
		if (missing(x0)) {
			x0 <- mean(e[1:2])
		}
		if (missing(y0)) {
			y0 <- mean(e[3:4])
		}
		angle <- angle[1]
		stopifnot(is.numeric(angle) && !is.nan(angle))
		x@ptr <- x@ptr$rotate(angle, x0[1], y0[1])
		messages(x, "spin")
	}
)
