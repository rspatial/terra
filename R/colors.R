

setMethod ("coltab" , "SpatRaster", 
	function(x) {
		hascols <- x@ptr$hasColors()
		if (any(hascols)) {
			d <- x@ptr$getColors()
			d <- lapply(d, .getSpatDF)
			d[!hascols] <- list(NULL)
		} else {
			d <- vector("list", length(hascols))
		}
		d
	}
)


setMethod ("coltab<-" , "SpatRaster", 
	function(x, layer=1, value) {

		x@ptr <- x@ptr$deepcopy()

		if (missing(value)) {
			value <- layer
			layer <- 1
		}
		layer <- layer[1]-1

		if (is.null(value)) {
			x@ptr$removeColors(layer)
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
		nms <- tolower(names(value))
		if (!("value" %in% nms)) {
			value <- cbind(values=(1:nrow(value))-1, value)
		}
		#value <- value[1:256,]
		if (ncol(value) == 4) {
			value <- cbind(value, alpha=255)
		}
		value[, 1] <- as.integer(value[, 1])
		for (i in 2:ncol(value)) {
			value[, i] <- as.integer(clamp(value[, i], 0, 255))
		} 
		value[is.na(value)] <- 255
		
		d <- .makeSpatDF(value)
		if (x@ptr$setColors(layer, d)) {
			return(x)
		} else {
			error("coltab<-", "cannot set these values")
		}
	}
)

