
setMethod("droplevels", signature(x="SpatRaster"),
	function(x, level=NULL, layer=1) {
		if (is.null(level)) {
			x@cpp <- x@cpp$droplevels()
			messages(x)
		} else {
			if (is.character(layer)) {
				layer <- match(layer, names(x))
				if (any(is.na(layer))) {
					error("droplevels", "invalid layer")
				}
			}
			x[[layer]][x[[layer]] %in%  level] <- NA
			x@cpp <- x@cpp$droplevels()
			messages(x)
		}
	}
)


setMethod("is.factor", signature(x="SpatRaster"),
	function(x) {
		x@cpp$hasCategories()
	}
)

setMethod("as.factor", signature(x="SpatRaster"),
	function(x) {
		x@cpp = x@cpp$makeCategorical(-1, spatOptions())
		messages(x)
		#if (!hasValues(x)) {
		#	error("as.factor", "x has no values")
		#}
		#x <- round(x)
		#u <- unique(x, TRUE)
		#for (i in 1:nlyr(x)) {
		#	set.cats(x, i, data.frame(ID=u[[i]], label=u[[i]], stringsAsFactors=FALSE))
		#}
		#x
	}
)


setMethod("levels", signature(x="SpatRaster"),
	function(x) {
		x <- x@cpp$getCategories()
		lapply(x, function(i) {
			d <- .getSpatDF(i$df)
			if (ncol(d) == 0) return("")
			d[, c(1, max(1, i$index+1))]
		})
	}
)


setMethod("levels<-", signature(x="SpatRaster"),
	function(x, value) {
		x@cpp <- x@cpp$deepcopy()
		if (is.null(value)) {
			x@cpp$removeCategories(-1)
			return(messages(x, "levels<-"))
		} else if (inherits(value, "list")) {
			for (i in 1:length(value)) {
				set.cats(x, i, value[[i]])
			}
		} else {
			set.cats(x, 1, value)
		}
		x
	}
)



setMethod ("set.cats" , "SpatRaster",
	function(x, layer=1, value, active=1) {

		if (missing(value)) {
			error("set.cats", "value cannot be missing")
			#return(invisible(x@cpp$setCatIndex(layer-1, index)))
		}

		if (is.character(layer)) {
			layer <- match(layer, names(x))
			if (any(is.na(layer))) {
				error("set.cats", "invalid layer")
			}
		}
		layer <- round(layer)

		if (length(layer) > 1) {
			if (!is.list(value)) {
				error("set.cats", "value should be a list")
			}
			if (length(layer) != length(value)) {
				error("set.cats", "length(value) != length(value)")
			}
			index <- rep_len(active, nlyr(x))
			for (i in 1:length(layer)) {
				ok <- set.cats(x, layer[i], value[[i]], index[i])
				x <- messages(x, "set.cats")
			}
			return(invisible(ok))
		} 

		if (layer < 1) {
			if (!is.list(value)) {
				error("set.cats", "value should be a list")
			}
			if (length(value) != nlyr(x)) {
				error("set.cats", "length(value) != nlyr(x)")
			}
			index <- rep_len(active, nlyr(x))
			for (i in 1:length(value)) {
				ok <- set.cats(x, i, value[[i]], index[i])
				x <- messages(x, "set.cats")
			}
			return(invisible(ok))
		}

		layer <- layer[1]
		if (is.character(layer)) {
			i <- match(layer, names(x))[1]
			if (length(i) == 0) {
				error("set.cats", layer, " is not in names(x)")
			}
			layer <- i
		} else {
			stopifnot(layer > 0 && layer <= nlyr(x))
		}

		if (inherits(value, "list")) {
			value <- value[[1]]
		}
		if (is.null(value)) {
			x@cpp$removeCategories(layer-1)
			messages(x, "set.cats")
			return(invisible(TRUE))
		}

		setname <- FALSE
		if (!is.data.frame(value)) {
			if (is.vector(value) || is.factor(value)) {
				if ((length(value) == 1) && value[1] == "") {
					return(invisible(""))
				}
				warn("set.cats", "setting categories like this is deprecated; use a two-column data.frame instead")
				value <- data.frame(value=0:(length(value)-1), category=value, stringsAsFactors=FALSE)
			} else {
				error("set.cats", "value should be a data.frame")
			}
		} else {
			setname <- TRUE
			if (ncol(value) == 1) {
				error("set.cats", "value should have at least two columns")
			} else {
				if (!is.numeric(value[[1]])) {
					error("set.cats", "the first column of 'value' must be numeric")
				}
				value[,1] <- round(value[[1]])
				if (length(unique(value[[1]])) != nrow(value)) {
					error("set.cats", "duplicate values (IDs) supplied")
				}
			}
		}
		value[[1]] <- as.integer(value[[1]])
		for (i in seq_along(value)) {
			if (is.factor(value[[i]])) {
				value[[i]] <- as.character(value[[i]])
			}
		}

		index <- max(1, min(ncol(value), active))
		if (setname) {
			nms <- names(x)
			cn <- colnames(value)[index+1]
			if (!(tolower(cn) %in% c("histogram", "count", "red", "green", "blue", "alpha", "opacity", "r", "g", "b", "a"))) {
				nms[layer] <- cn
				if (! x@cpp$setNames(nms, FALSE)) {
					error("names<-", "cannot set name")
				}
			}
		}
		if (any(is.na(value[[1]]))) {
			error("set.cats", "you cannot associate a category with NA")
		}
		if (any(table(value[[1]]) > 1)) {
			error("set.cats", "you cannot have duplicate IDs")
		}

		value <- .makeSpatDF(value)
		ok <- x@cpp$setCategories(layer-1, value, index)
		x <- messages(x, "set.cats")
		invisible(ok)
	}
)



setMethod ("categories" , "SpatRaster",
	function(x, layer=1, value, active=1, ...) {
		#... to accept but ignore old argument "index"
		x@cpp <- x@cpp$deepcopy()
		set.cats(x, layer, value, active)
		x
	}
)


setMethod ("activeCat" , "SpatRaster",
	function(x, layer=1) {
		layer <- layer[1]
		if (is.character(layer)) {
			layer = which(layer == names(x))[1]
			if (is.na(layer)) {
				error("activeCat", "invalid layer name")
			}
		}
		if (layer < 1) {
			sapply(1:nlyr(x), function(i) x@cpp$getCatIndex(i-1))
		} else {
			if (!is.factor(x)[layer]) {
				return(NA)
			}
			x@cpp$getCatIndex(layer-1)
		}
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
		x <- deepcopy(x)
		if (!x@cpp$setCatIndex(layer-1, value)) {
			error("activeCat", "invalid category index")
		}
		x
	}
)

setMethod("cats" , "SpatRaster",
	function(x, layer) {
		if (!missing(layer)) {
			x <- subset(x, layer, NSE=FALSE)
		}
		cats <- x@cpp$getCategories()
		lapply(1:nlyr(x), function(i) {
			if (cats[[i]]$df$nrow == 0) {
				return(NULL)
			}
			.getSpatDF(cats[[i]]$df)
		})
	}
)


# superseded by levels(x)[[layer]]
..active_cats <- function(x, layer) {
	ff <- is.factor(x)
	if (!any(ff)) {
		return (lapply(ff, function(i) NULL))
	}
	cats <- x@cpp$getCategories()
	x <- lapply(1:length(cats), function(i) {
		if (cats[[i]]$df$nrow == 0) return(NULL)
		r <- .getSpatDF(cats[[i]]$df)
		a <- activeCat(x, i)
		if (a < 0) return(NULL)
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
		if (!any(is.factor(x))) {
			x <- deepcopy(x)
			x@cpp$setValueType(0)
			return(x)
		}
		if (nlyr(x) > 1) {
			x <- lapply(1:nlyr(x), function(i) as.numeric(x[[i]], index=index))
			x <- rast(x)
			if (filename != "") {
				x <- writeRaster(x, filename, ...)
			}
			return(x)
		}
		g <- cats(x)[[1]]
		if (!is.null(index)) {
			if (is.character(index)) {
				index <- match(index, colnames(g))
				if (is.na(index)) {
					error("as.numeric", "index is not category name")				
				}
				if (index == 1) {
					levels(x) <- NULL
					x@cpp$setValueType(0)
					if (filename != "") {
						x <- writeRaster(x, filename, ...)
					}
					return(x)
				}		
			} else {
				index <- round(index[1])
				if (!((index >= 1) & (index < ncol(g)))) {
					error("as.numeric", "index out of range")
				}
				index <- index + 1
			}
		} else {
			index <- activeCat(x, 1)
			if (index <= 1) {
				levels(x) <- NULL
				x@cpp$setValueType(0)
				if (filename != "") {
					x <- writeRaster(x, filename, ...)
				}
				return(x)
			}
		}
		from <- g[,1]
		to <- g[,index]
		if (!is.numeric(to)) {
			suppressWarnings(toto <- as.numeric(to))
			if (sum(is.na(toto) > sum(is.na(to)))) {
				to <- as.integer(as.factor(to))
			} else {
				to <- toto
			}
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
			set.cats(x, 1, fact)
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
					z <- as.numeric(y, index=j-1)
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



setMethod("concats", "SpatRaster",
	function(x, y, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@cpp = x@cpp$combineCats(y@cpp, opt)
		messages(x, "concats")
	}
)

setMethod("addCats", "SpatRaster",
	function(x, value, merge=FALSE, layer=1) {
		if (!(is.factor(x)[layer])) {
			error("addCat", "the layer has no categories to add to")
		}
		cts <- cats(x)[[layer]]
		nact <- ncol(cts)
		if (merge) {
			if (ncol(value) < 2) {
				error("addCat", "'value' must have at least two columns when using 'merge=TRUE'")
			}
			cts <- merge(cts, value, by=1, all.x=TRUE)
			cts <- cts[order(cts[,1]), ]
		} else {
			if (nrow(cts) != nrow(value)) {
				error("addCat", "the number of categories does not match")
			}
			cts <- cbind(cts, value)
		}
		categories(x, layer=layer, cts, active=nact)
	}
)

