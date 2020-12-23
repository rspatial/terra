

setMethod("intersect", signature(x="SpatVector", y="SpatVector"), 
	function(x, y) {
		x <- methods::as(x, "Spatial")
		y <- methods::as(y, "Spatial")
		r <- raster::intersect(x, y)
		vect(r)
	}
)

## cheating --- using raster/sp/rgeos for now
setMethod("crop", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, ...) {
		x <- methods::as(x, "Spatial")
		y <- methods::as(y, "Spatial")
		r <- raster::crop(x, y)
		vect(r)
	}
)

setMethod("crop", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y, ...) {
		y <- as.polygons(y)
		crop(x, y, ...)
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



setMethod("disaggregate", signature(x="SpatVector"), 
	function(x, ...) {
		x@ptr <- x@ptr$disaggregate()
		messages(x, "disaggregate")
	}
)




