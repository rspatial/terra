
setMethod("buffer", signature(x="SpatVector"), 
	function(x, width, quadsegs=10, capstyle="round", ...) {
		if (geomtype(x) == "points") {
			x@ptr <- x@ptr$buffer(width, quadsegs, 1)
		} else {
			if (isLonLat(x)) {
				warning("lon/lat data are treated as planar")
			}
			x@ptr <- x@ptr$buffer2(width, quadsegs, 1)		
		}
		show_messages(x, "buffer")
	}
)


## cheating --- using raster/rgeos for now

setMethod("crop", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, ...) {
		x <- methods::as(x, "Spatial")
		y <- methods::as(y, "Spatial")
		r <- raster::crop(x, y)
		vect(r)
	}
)


setMethod("aggregate", signature(x="SpatVector"),
	function(x, by=NULL, sums=NULL, dissolve=TRUE, vars=NULL, ...) {
		#gt <- geomtype(x)
		x <- methods::as(x, "Spatial")
		r <- aggregate(x, by=by, sums=sums, dissolve=dissolve, vars=vars, ...)
		vect(r)
	}
)


setMethod("intersect", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x <- methods::as(x, "Spatial")
		y <- methods::as(y, "Spatial")
		r <- raster::intersect(x, y)
		vect(r)
	}
)

