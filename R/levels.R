

setMethod("is.factor", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasCategories()
	}
)

setMethod("as.factor", signature(x="SpatRaster"), 
	function(x) {
		x <- round(x)
		u <- unique(x)
		levels(x) <- cbind(u, u)
		x
	}
)


setMethod("levels", signature(x="SpatRaster"), 
	function(x) {
		x <- x@ptr$getCategories()
		lapply(x, function(i) {
			d <- .getSpatDF(i$df)
			if (ncol(d) == 0) return("")
			d[, c(1, max(1, i$index+1))]
		})
	}
)


setMethod("levels<-", signature(x="SpatRaster"), 
	function(x, value) {
		x@ptr <- x@ptr$deepcopy()
		if (is.null(value)) {
			x@ptr$removeCategories(0)
			return(messages(x, "levels<-"))
		} else if (inherits(value, "list")) {
			for (i in 1:length(value)) {
				set.cats(x, i, value[[i]])
			}
		} else {
			set.cats(x, 1, value, 2)
		}
		x
	}
)




setMethod ("set.cats" , "SpatRaster", 
	function(x, layer=1, value, index) {
		layer <- layer[1]
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
			messages(x, "set.cats")
			return(invisible(TRUE))
		}

		if (inherits(value, "list")) {
			value <- value[[1]]
		}
		setname <- FALSE
		if (!is.data.frame(value)) {
			if (is.vector(value) || is.factor(value)) {
				if ((length(value) == 1) && value[1] == "") {
					return(invisible(""))
				}
				warning("set.cats", "setting categories like this is deprecated; use a two-column data.frame instead")
				value <- data.frame(value=0:(length(value)-1), category=value)
			} else {
				error("set.cats", "value should be a data.frame or a vector")
			}
		} else {
			setname <- TRUE
			if (ncol(value) == 1) {
				value <- data.frame(value=1:nrow(value), value)
			} else {
				value[,1] <- round(value[,1])
				if (length(unique(value[,1])) != nrow(value)) {
					error("set.cats", "duplicate values (IDs) supplied")
				}
			}
		}

		index <- max(1, min(ncol(value), index))
#		if (is.data.frame(value)) {
		if (setname) {
			nms <- names(x)
			nms[layer] <-  colnames(value)[index]
			if (! x@ptr$setNames(nms, FALSE)) {
				error("names<-", "cannot set name")
			}
		}
		if (any(is.na(value[,1]))) {
			error("set.cats", "you cannot associate a category with NA")
		}

		value <- .makeSpatDF(value)
		ok <- x@ptr$setCategories(layer-1, value, index-1)
#		} else {
#			value <- as.character(value)
#			x@ptr$setLabels(layer-1, value)
#		}
		x <- messages(x, "set.cats")
		invisible(ok)
	}
)


setMethod ("categories" , "SpatRaster", 
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

		x@ptr <- x@ptr$deepcopy()
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
			return(messages(x, "set.cats"))
		}


		if (inherits(value, "list")) {
			value <- value[[1]]
		}
		setname <- FALSE
		if (!is.data.frame(value)) {
			if (is.vector(value) || is.factor(value)) {
				if ((length(value) == 1) && value[1] == "") {
					return(invisible(""))
				}
				value <- data.frame(value=0:(length(value)-1), category=value)
			} else {
				error("set.cats", "value should be a data.frame or a vector")
			}
		} else {
			setname <- TRUE
			if (ncol(value) == 1) {
				value <- data.frame(value=1:nrow(value), value)
			} else {
				value[,1] <- round(value[,1])
				if (length(unique(value[,1])) != nrow(value)) {
					error("setCats", "duplicate ID values supplied")
				}
			}
		}
		#	v <- data.frame(value=0:maxv)
		#	value <- merge(v, value, by=1, all.x=TRUE)
		#}

		index <- max(1, min(ncol(value), index))
#		if (is.data.frame(value)) {
		if (setname) {
			nms <- names(x)
			nms[layer] <-  colnames(value)[index]
			if (! x@ptr$setNames(nms, FALSE)) {
				error("names<-", "cannot set name")
			}
		}

		if (any(is.na(value[,1]))) {
			error("categories", "you cannot associate a category with NA")
		}
		value <- .makeSpatDF(value)
		ok <- x@ptr$setCategories(layer-1, value, index-1)
#		} else {
#			value <- as.character(value)
#			x@ptr$setLabels(layer-1, value)
#		}
		messages(x, "categories")
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
	function(x, layer, active=FALSE) {
		if (!missing(layer)) {
			x <- subset(x, layer)
		}
		cats <- x@ptr$getCategories()
		lapply(1:nlyr(x), function(i) {
			if (cats[[i]]$df$nrow == 0) {
				return(NULL)
			}
			y <- .getSpatDF(cats[[i]]$df)
			if (active) {
				y <- y[, c(1, activeCat(x[[i]], i) + 1)]
			}
			y
		})
	}
)



active_cats <- function(x, layer) {
	ff <- is.factor(x)
	if (!any(ff)) {
		return (lapply(ff, function(i) NULL))
	}
	cats <- x@ptr$getCategories()
	x <- lapply(1:length(cats), function(i) {
		if (cats[[i]]$df$nrow == 0) return(NULL)
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
			index <- set.cats(x, 1)
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
			index <- set.cats(x, 1)
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
			set.cats(x, 1, fact, 2)
		}
		x
}



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


