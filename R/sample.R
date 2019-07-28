

setMethod("sampleRaster", signature(x="SpatRaster", size="numeric"), 
	function(x, size, type="regular", ...) { 
		size <- max(1, min(size(x), size))
		#if (type == "regular") {
			x@ptr <- x@ptr$sampleRegular(size)
		#}
		show_messages(x, "sampleRegular")		
	}
)


setMethod("sampleCells", signature(x="SpatRaster", size="numeric"), 
	function(x, size, type="regular", xy=FALSE, vect=FALSE, ...) { 
		type <- tolower(type)
		size <- max(1, min(size(x), size))
		if (vect) xy=TRUE
		if (type == "regular") {
			x@ptr <- x@ptr$sampleRegular(size)
			x <- show_messages(x, "sampleRegular")	
			if (xy) {
				pts <- xyFromCell(x, 1:ncell(x))
				if (vect) {
					v <- vect(pts, "points", data.frame(values(x)), crs=crs(x))
					return(v)
				} else {
					pts <- cbind(pts, values(x))
					colnames(pts)[1:2] <- c("x", "y")
					return(pts)
				}
			} else {
				return(values(x))
			}
		} else {
			stop("type unknown or not yet implemented")
		}
	}
)
