
setMethod ("has.colors" , "SpatRaster",
	function(x) {
		x@pntr$hasColors()
	}
)



setMethod ("coltab" , "SpatRaster",
	function(x) {
		hascols <- x@pntr$hasColors()
		if (any(hascols)) {
			d <- x@pntr$getColors()
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
		
		x@pntr <- x@pntr$deepcopy()
		if (inherits(value, "list")) {
			for (i in seq_along(value)) {
				layer <- i-1
				if (is.null(value[[i]])) {
					x@pntr$removeColors(layer)
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
				if (!x@pntr$setColors(layer, d)) {
					messages(x, "cols<-")
				}
			}
		} else {
			if (inherits(layer, "character")) {
				layer <- match(layer, names(x))
				if (is.na(layer)) {
					error("coltab", "not a valid layer name")
				}
			}
			layer <- layer[1]-1
			if (is.null(value)) {
				x@pntr$removeColors(layer)
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
			} #else {
			#	nms <- tolower(names(value))
			#	if (!(grepl("value", nms))) {
			#		value <- cbind(values=(1:nrow(value))-1, value)
			#	}
			#	#value <- value[1:256,]
			#	if (ncol(value) == 4) {
			#		value <- cbind(value, alpha=255)
			#	}
			#}
			value[, 1] <- as.integer(value[, 1])
			for (i in 2:ncol(value)) {
				value[, i] <- as.integer(clamp(value[, i], 0, 255))
			}
			value[is.na(value)] <- 255
			d <- .makeSpatDF(value)
			if (!x@pntr$setColors(layer, d)) {
				messages(x, "cols<-")
			}
		}
		x
	}
)




make.RGB <- function(x, col=grDevices::rainbow(25), breaks=NULL, alpha=FALSE, colNA="white", zlim=NULL, zlimcol=NULL, ext=NULL, filename="", ...) { 

	getCols <- function(x, col, breaks=NULL, r=NULL, colNA=NA) {
		if (!is.null(breaks)) {
			breaks <- sort(breaks)
			x <- as.numeric(cut(x, breaks, include.lowest=TRUE))
			
		} else {
			x <- (x - r[1])/ (r[2] - r[1])
			x <- round(x * (length(col)-1) + 1)
		}
		x <- col[x]
		if (!is.na(colNA)) {
			x[is.na(x)] <- grDevices::rgb(t(grDevices::col2rgb(colNA)), maxColorValue=255)
		}
		x
	}

	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	
	if (alpha) {
		out <- rast(x, nl=4)
		RGB(out) <- 1:4
	} else {
		out <- rast(x, nl=3)
		RGB(out) <- 1:3
	}
	names(out) <- c('red', 'green', 'blue', 'alpha')[1:nlyr(out)]

	
	r <- minmax(x, TRUE)[,1]
	if (is.null(breaks)) {
		zrange <- range(r, zlim, na.rm=TRUE)
	} else {
		zrange <- range(r, zlim, breaks, na.rm=TRUE)
	}
	if (zrange[1] == zrange[2]) {
		zrange[1] <- zrange[1] - 0.001
		zrange[2] <- zrange[2] + 0.001
	}

	tr <- writeStart(out, filename=filename, ...)
		
	for (i in 1:tr$n) {
		v <- values(x, row=tr$row[i], nrows=tr$nrows[i])
		
		if (!is.null(zlim)) {
			if (!is.null(zlimcol)) {
				v[v < zlim[1]] <- zlim[1]
				v[v > zlim[2]] <- zlim[2]
			} else { #if (is.na(zlimcol)) {
				v[v < zlim[1] | v > zlim[2]] <- NA
			} 
		}
		v <- getCols(v, col, breaks, zrange, colNA)
		v <- grDevices::col2rgb(as.vector(v), alpha=alpha)
		writeValues(out, t(v), tr$row[i], tr$nrows[i])
	}
	writeStop(out)
}
