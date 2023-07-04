
setMethod ("has.colors" , "SpatRaster",
	function(x) {
		x@pnt$hasColors()
	}
)



setMethod ("coltab" , "SpatRaster",
	function(x) {
		hascols <- x@pnt$hasColors()
		if (any(hascols)) {
			d <- x@pnt$getColors()
			d <- lapply(d, .getSpatDF)
			d[!hascols] <- list(NULL)
		} else {
			d <- vector("list", length(hascols))
		}
		d
	}
)


setMethod ("coltab<-" , "SpatRaster",
	function(x, ..., layer=1, value) {

		x@pnt <- x@pnt$deepcopy()
		if (inherits(value, "list")) {
			for (i in seq_along(value)) {
				layer <- i
				if (is.null(value[[i]])) {
					x@pnt$removeColors(layer)
				}
				if (inherits(value[[i]], "character")) {
					value[[i]] <- data.frame(t(grDevices::col2rgb(value[[i]], alpha=TRUE)), stringsAsFactors=FALSE)
				} else if (inherits(value[[i]], "matrix")) {
					value[[i]] <- data.frame(value[[i]])
				}
				if (!inherits(value[[i]], "data.frame")) {
					error("coltab<-", "cannot process these color values")
				}
				if (ncol(value[[i]]) == 2) {
					value[[i]] <- data.frame(values=value[[i]][,1], t(grDevices::col2rgb(value[[i]][,2], alpha=TRUE)), stringsAsFactors=FALSE)
				} 
				value[[i]][, 1] <- as.integer(value[[i]][, 1])
				for (j in 2:ncol(value[[i]])) {
					value[[i]][, j] <- as.integer(clamp(value[[i]][, j], 0, 255))
				}
				value[[i]][is.na(value[[i]])] <- 255

				d <- .makeSpatDF(value[[i]])
				if (!x@pnt$setColors(layer, d)) {
					error("coltab<-", "cannot set these values")
				}
			}		
		} else {
			layer <- layer[1]-1
			if (is.null(value)) {
				x@pnt$removeColors(layer)
				return(x)
			}

			if (inherits(value, "character")) {
				value <- data.frame(t(grDevices::col2rgb(value, alpha=TRUE)), stringsAsFactors=FALSE)
			} else if (inherits(value, "matrix")) {
				value <- data.frame(value)
			}
			if (!inherits(value, "data.frame")) {
				error("coltab<-", "cannot process these color values")
			}
			if (ncol(value) == 2) {
				value <- data.frame(values=value[,1], t(grDevices::col2rgb(value[,2], alpha=TRUE)), stringsAsFactors=FALSE)
			} else {
			#	nms <- tolower(names(value))
			#	if (!(grepl("value", nms))) {
			#		value <- cbind(values=(1:nrow(value))-1, value)
			#	}
			#	#value <- value[1:256,]
			#	if (ncol(value) == 4) {
			#		value <- cbind(value, alpha=255)
			#	}
			}
			value[, 1] <- as.integer(value[, 1])
			for (i in 2:ncol(value)) {
				value[, i] <- as.integer(clamp(value[, i], 0, 255))
			}
			value[is.na(value)] <- 255
			d <- .makeSpatDF(value)
			if (!x@pnt$setColors(layer, d)) {
				messages("cols<-", x)
			}
			x
		}
	}
)

