

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
				if (length(value == 1) && value[1] == "") {
					return(invisible(""))
				}
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


setMethod ("activeCat" , "SpatRaster", 
	function(x, layer=1) {
		layer = layer[1]
		if (is.character(layer)) {
			layer = which(layer == names(x))[1]
			if (is.na(layer)) {
				error("activeCat", "invalid layer name")
			}
		}
		if (!is.factor(x)[layer]) {
			return(NA)
		}
		x@ptr$getCatIndex(layer-1)
	}
)

setMethod("activeCat<-" , "SpatRaster", 
	function(x, layer=1, value) {
		if (missing(value)) {
			value <- layer[1]
			layer <- 1
		} else {
			layer <- layer[1]
		}
		if ((layer < 1) | (layer > nlyr(x))) {
			error("activeCat", "invalid layer")		
		}
		if (!is.factor(x)[layer]) {
			error("activeCat", "layer is not categorical")				
		}
		if (is.character(value)) {
			g <- cats(x)[[layer]]
			value <- which(value == names(g))[1] - 1
			if (is.na(value)) {
				error("activeCat", "invalid category name")
			}
		}
		if (!x@ptr$setCatIndex(layer-1, value)) {
			error("activeCat", "invalid category index")
		} 
		x
	}
)

setMethod("cats" , "SpatRaster", 
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



active_cats <- function(x, layer) {
	ff <- is.factor(x)
	if (!any(ff)) {
		return (lapply(ff, function(i)NULL))
	}
	cats <- x@ptr$getCategories()
	x <- lapply(1:length(cats), function(i) {
		if (cats[[1]]$df$nrow == 0) return(NULL)
		r <- .getSpatDF(cats[[i]]$df)
		a <- activeCat(x, i)
		r[, c(1, a+1)]
	})
	
	if (!missing(layer)) {
		x[[layer]]
	} else {
		x
	}
}





setMethod ("as.numeric", "SpatRaster", 
	function(x, index=NULL, filename="", ...) {
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
			to <- as.integer(as.factor(to))
		}
		m <- cbind(from, to)
		m <- m[!is.na(m[,1]), ,drop=FALSE]
		classify(x, m, names=names(g)[index], filename, ...)
	}
)

	

catLayer <- function(x, index, ...) {
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
		toc <- g[,index]
		
		addFact <- FALSE
		if (!is.numeric(toc)) {
			addFact <- TRUE
			ton <- as.integer(as.factor(toc))
		} else {
			ton <- toc
		}
		m <- cbind(from, ton)
		m <- m[!is.na(m[,1]), ,drop=FALSE]
		x <- classify(x, m, names=names(g)[index], ...)
		if (addFact) {
			fact <- unique(data.frame(ton, toc))
			names(fact) <- c("ID", names(g)[index])
			fact <- fact[order(fact[,1]), ]
			setCats(x, 1, fact, 2)
		}
		x
	}
)
	
	
setMethod("catalyze", "SpatRaster", 
	function(x, filename="", ...) {
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
		out <- rast(out)
		if (filename!="") {
			out <- writeRaster(out, filename, ...)
		}
		out
	}
)


