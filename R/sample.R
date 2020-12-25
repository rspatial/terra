
.get_seed <- function() {
  sample.int(.Machine$integer.max, 1)
}


.sampleCells <- function(x, size, method, replace) {
	if (method == "random") {
		cells <- sample(ncell(x), size, replace=replace)
	} else { # regular 
		f <- sqrt(size / ncell(x))
		nr <- ceiling(nrow(x) * f)
		nc <- ceiling(ncol(x) * f);
		xstep <- ncol(x) / nc
		ystep <- nrow(x) / nr
		xsamp <- seq(0.5*xstep, ncol(x), xstep)
		ysamp <- seq(0.5*ystep, nrow(x), ystep)
		xy <- expand.grid(round(ysamp), round(xsamp))
		cells <- cellFromRowCol(x, xy[,1], xy[,2]) 
	}
	return(cells)
}

setMethod("spatSample", signature(x="SpatRaster"), 
	function(x, size, method="regular", replace=FALSE, as.raster=FALSE, cells=FALSE, ...) {

		if (cells) {
			return(.sampleCells(x, size, method, replace))
		}
		if (!hasValues(x) & !as.raster) {
			error("spatSample", "SpatRaster has no values. Use cells=TRUE or as.raster=TRUE")
		}
		method <- tolower(method)
		stopifnot(method %in% c("random", "regular"))
		size <- round(size)
		stopifnot(size > 0)
		size <- min(ncell(x), size)

		if (method == "regular") {
			if (as.raster) {
				x@ptr <- x@ptr$sampleRegularRaster(size)
				x <- messages(x, "spatSample")
				return(x);
			} else {
				v <- x@ptr$sampleRegularValues(size)
			}
		} else {
			seed <- .get_seed()
			if (as.raster) {
				x@ptr <- x@ptr$sampleRandomRaster(size, replace, seed)
				x <- messages(x, "spatSample")
				return(x);
			} else {
				v <- x@ptr$sampleRandomValues(size, replace, seed)
			}
		}
		# values
		x <- messages(x, "spatSample")
		if (length(v) > 0) {
			v <- do.call(cbind, v)
			colnames(v) <- names(x)
		}
		return(v)
	}
)



#setMethod("spatSample", signature(x="SpatRaster", size="numeric"), 
#	function(x, size, ...) { 
#		size <- max(1, min(size(x), size))
#		x@ptr <- x@ptr$spatSample(size)
#		messages(x, "spatSample")
#	}
#)


# setMethod("sampleCells", signature(x="SpatRaster", size="numeric"), 
	# function(x, size, type="regular", xy=FALSE, vect=FALSE, ...) { 
		# type <- tolower(type)
		# size <- max(1, min(size(x), size))
		# if (vect) xy=TRUE
		# if (type == "regular") {
			# x@ptr <- x@ptr$spatSample(size)
			# x <- messages(x, "spatSample")
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
