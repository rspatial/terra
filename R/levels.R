

setMethod("is.factor", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasCategories()
	}
)


setMethod("levels", signature(x="SpatRaster"), 
	function(x) {
		x <- x@ptr$getCategories()
		lapply(x, function(i) {
			d <- .getSpatDF(i$df)
			if (ncol(d) == 0) return("")
			d[,max(1, i$index+1)]
		})
	}
)


setMethod("levels<-", signature(x="SpatRaster"), 
	function(x, value) {
		if (is.null(value)) {
			x@ptr$removeCategories(0)
			return(messages(x, "levels<-"))
		} else if (inherits(value, "list")) {
			for (i in 1:length(value)) {
				setCats(x, i, value[[i]])
			}
		} else {
			setCats(x, 1, value, 2)		
		}
		x
	}
)


setMethod ("setCats" , "SpatRaster", 
	function(x, layer=1, value, index) {
		layer = layer[1]
		if (is.character(layer)) {
			i <- match(layer, names(x))[1]
			if (length(i) == 0) {
				error("setLevels", layer, " is not in names(x)")
			}
			layer <- i
		} else {
			stopifnot(layer > 0 && layer <= nlyr(x))
		}

		if (missing(value)) {
			if (missing(index)) {
				return(x@ptr$getCatIndex(layer-1) + 1)
			} else {
				return(invisible(x@ptr$setCatIndex(layer-1, index)))
			}
		} 
		if (missing(index)) {
			index <- 2
		}
		if (is.null(value)) {
			x@ptr$removeCategories(layer-1)
			return(messages(x, "setCats"))
		}


		if (inherits(value, "list")) {
			value <- value[[1]]		
		}
		setname <- FALSE
		if (!is.data.frame(value)) {
			if (is.vector(value) || is.factor(value)) {
				value <- data.frame(ID=0:(length(value)-1), category=value)
			} else {
				error("setCats", "value should be a data.frame or a vector")
			}
		} else {
			setname <- TRUE
			if (nrow(value) > 256) {
				error("setCats", "you can set no more than 256 categories")
			}
			if (ncol(value) == 1) {
				value <- data.frame(ID=1:nrow(value), value)
			} else {
				value[,1] <- round(value[,1])
				if (length(unique(value[,1])) != nrow(value)) {
					error("setCats", "duplicate IDs supplied")
				}
				r <- range(value[,1])
				if (r[1] < 0 || r[2] > 255) {
					error("seCats", "ID values must be between 0 and 255")
				}
			}
		}
		maxv <- max(value[,1])
		v <- data.frame(ID=0:maxv)
		value <- merge(v, value, by=1, all.x=TRUE)
			
		
		index <- max(1, min(ncol(value), index))
#		if (is.data.frame(value)) {
		if (setname) {
			names(x)[layer] <- colnames(value)[index]
		}

		value <- .makeSpatDF(value)
		ok <- x@ptr$setCategories(layer-1, value, index-1)
#		} else {
#			value <- as.character(value)
#			x@ptr$setLabels(layer-1, value)
#		}
		x <- messages(x, "setCats")
		invisible(ok)
	}
)


setMethod ("cats" , "SpatRaster", 
	function(x, layer) {
		x <- x@ptr$getCategories()
		x <- lapply(x, function(i) {
			if (is.null(i)) return( NULL)
			.getSpatDF(i$df)
		})
		if (!missing(layer)) {
			x[[layer]]
		} else {
			x
		}
	}
)




setMethod ("as.numeric" , "SpatRaster", 
	function(x, index=NULL, ...) {
		stopifnot(nlyr(x) == 1)
		if (!is.factor(x)) return(x)
		g <- cats(x)[[1]]
		if (!is.null(index)) {
			if (!((index > 1) & (index <= ncol(g)))) {
				error("as.numeric", "invalid index")
			}
		} else {
			index <- setCats(x, 1)
		}
		from <- g[,1]
		to <- g[,index]
		if (!is.numeric(to)) {
			to <- 1:length(to)
		}
		x <- classify(x, cbind(from, to), names=names(g)[index], ...)
		messages(x)
	}
)


	
	
.as.layers <- function(x) {
	g <- cats(x)
	out <- list()
	for (i in 1:nlyr(x)) {
		y <- x[[i]]
		gg <- g[[i]]
		if (nrow(gg) > 0) {
			for (j in 2:ncol(gg)) {
				z <- as.numeric(y, index=j)
				out <- c(out, z)
			}
		} else {
			out <- c(out, y)
		}
	}
	rast(out)
}

