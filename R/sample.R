
.get_seed <- function() {
  sample.int(.Machine$integer.max, 1)
}

setMethod("spatSample", signature(x="SpatRaster"), 
	function(x, size, method="regular", replace=FALSE, as.raster=FALSE, ...) {
		method <- tolower(method)
		stopifnot(method %in% c("random", "regular"))
		size <- round(size)
		stopifnot(size > 0)
		size <- min(ncell(x), size)
		
		if (method == "regular") {
			if (as.raster) {
				x@ptr <- x@ptr$sampleRegular(size)
				x <- show_messages(x, "spatSample")		
				return(x);
			} else {
				x@ptr <- x@ptr$sampleRegular(size)
				x <- show_messages(x, "spatSample")		
				return(values(x))
			}
		} else {
			seed <- .get_seed()
			r <- x@ptr$sampleRandom(size, replace, seed)
			show_messages(x, "spatSample")		
			return(r)
		}
	}
)



#setMethod("spatSample", signature(x="SpatRaster", size="numeric"), 
#	function(x, size, ...) { 
#		size <- max(1, min(size(x), size))
#		x@ptr <- x@ptr$spatSample(size)
#		show_messages(x, "spatSample")		
#	}
#)


# setMethod("sampleCells", signature(x="SpatRaster", size="numeric"), 
	# function(x, size, type="regular", xy=FALSE, vect=FALSE, ...) { 
		# type <- tolower(type)
		# size <- max(1, min(size(x), size))
		# if (vect) xy=TRUE
		# if (type == "regular") {
			# x@ptr <- x@ptr$spatSample(size)
			# x <- show_messages(x, "spatSample")	
			# if (xy) {
				# pts <- xyFromCell(x, 1:ncell(x))
				# if (vect) {
					# v <- vect(pts, "points", data.frame(values(x)), crs=crs(x))
					# return(v)
				# } else {
					# pts <- cbind(pts, values(x))
					# colnames(pts)[1:2] <- c("x", "y")
					# return(pts)
				# }
			# } else {
				# return(values(x))
			# }
		# } else {
			# stop("type unknown or not yet implemented")
		# }
	# }
# )
