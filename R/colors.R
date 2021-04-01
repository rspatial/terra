

setMethod ("coltab" , "SpatRaster", 
	function(x) {
		hascols <- x@ptr$hasColors()
		if (any(hascols)) {
			d <- x@ptr$getColors()
			d <- lapply(d, .getSpatDF)
		} else {
			d <- vector("list", length(hascols))
		}
		d
	}
)


setMethod ("coltab<-" , "SpatRaster", 
	function(x, layer=1, value) {
		stopifnot(hasValues(x))
		if (missing(value)) {
			value <- layer
			layer <- 1
		}
		layer <- layer[1]-1
		
		if (is.null(value)) {
			x@ptr$removeColors(layer-1)
			return(x)
		}

		if (inherits(value, "list")) {
			value <- value[[1]]
		}
		if (inherits(value, "character")) {
			value <- t(grDevices::col2rgb(value, alpha=TRUE))
		}
		if (inherits(value, "character")) {
			value <- data.frame(t(grDevices::col2rgb(value, alpha=TRUE)))
		} else if (inherits(value, "matrix")) {
			value <- data.frame(value)
		}

		stopifnot(inherits(value, "data.frame"))

		value <- value[1:256,]
		value[is.na(value)] <- 255
		value <- data.frame(sapply(value, function(i) as.integer(clamp(i, 0, 255))))
		if (ncol(value) == 3) {
			value <- cbind(value, alpha=255)
		}

		d <- .makeSpatDF(value)
		if (x@ptr$setColors(layer, d)) {
			return(x)
		} else {
			error("coltab<-", "cannot set these values")
		}
	}
)

