

setMethod("is.factor", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasCategories()
	}
)


setMethod("levels", signature(x="SpatRaster"), 
	function(x) {
		x <- x@ptr$getCategories()
		lapply(x, function(i) {
			if (is.null(i)) return( NULL)
			d <- .getSpatDF(i$df)
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
				setCats(x, i, value[[i]], 1)
			}
		} else {
			setCats(x, 1, value, 1)		
		}
		x
	}
)


setMethod ("setCats" , "SpatRaster", 
	function(x, layer=1, value, index=1) {
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
		if (!is.data.frame(value)) {
			if (is.list(value)) {
				value <- value[[1]]
			}
		}

		if (is.null(value)) {
			x@ptr$removeCategories(layer-1)
			return(messages(x, "setLevels"))
		}
		if (is.data.frame(value)) {
			value <- .getSpatDF(value)
			x@ptr$setCategories(layer-1, value, index-1)
		} else {
			value <- as.character(value)
			x@ptr$setLabels(layer-1, value)
		}
		messages(x, "setLevels")
	}
)


setMethod ("cats" , "SpatRaster", 
	function(x) {
		x <- x@ptr$getCategories()
		lapply(x, function(i) {
			if (is.null(i)) return( NULL)
			.getSpatDF(i$df)
		})
	}
)

catvar <- function(x, layer=1) {
	stopifnot(layer > 0)
	x@ptr$getCatIndex(layer-1);
}

setcatvar <- function(x, layer=1, index=1) {
	stopifnot(layer > 0)
	stopifnot(index > 0)
	x@ptr$setCatIndex(layer-1, index-1);
}


