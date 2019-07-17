
setMethod("buffer", signature(x="SpatVector"), 
	function(x, width, quadsegs=10, capstyle="round", ...) {
		if (geomtype(p) == "points") {
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

