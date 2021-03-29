
setMethod ("rats" , "SpatRaster", 
	function(x) {
		att <- x@ptr$hasAttributes()
		if (any(att)) {
			d <- x@ptr$getAttributes()
			d <- lapply(d, .getSpatDF)
		} else {
			d <- vector("list", length(att))
		}
		d
	}
)


setMethod ("rats<-" , "SpatRaster", 
	function(x, layer=1, value) {
		stopifnot(is.data.frame(value))
		if (missing(value)) {
			value <- layer
			layer <- 1
		}
		if (is.character(layer)) {
			i <- match(layer, names(x))[1]
			if (length(i) == 0) {
				error("rats<-", layer, " is not in names(x)")
			}
			layer <- i
		} else {
			stopifnot(layer > 0 && layer <= nlyr(x))
		}
		rat <- .makeSpatDF(value)
		x@ptr$setAttributes(layer-1, value)
		x
	}
)



setMethod("as.factor", signature(x="SpatRaster"), 
	function(x) {
		stopifnot(hasValues(x))
		opt <- spatOptions()
		x@ptr <- x@ptr$makeCategorical(0, opt)
		messages(x, "as.factor")
	}
)

setMethod("is.factor", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasCategories()
	}
)


setMethod("levels", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$getCategories()
	}
)

setMethod("cats", signature(x="SpatRaster"), 
	function(x) {
		levels(x)
	}
)


setMethod("levels<-", signature(x="SpatRaster"), 
	function(x, value) {
		cats(x, 1) <- value
		x
	}
)


setMethod ("cats<-" , "SpatRaster", 
	function(x, layer=1, value) {
		stopifnot(hasValues(x))
		stopifnot(is.data.frame(value))
		if (missing(value)) {
			value <- layer
			layer <- 1
		}
		if (is.character(layer)) {
			i <- match(layer, names(x))[1]
			if (length(i) == 0) {
				error("cats<-", layer, " is not in names(x)")
			}
			layer <- i
		} else {
			stopifnot(layer > 0 && layer <= nlyr(x))
		}
		if (is.null(value) || is.na(value[[1]][1])) {
			x@ptr$removeCategories(layer-1)
			return(messages(x, "levels<-"))
		}
		if (!is.factor(x[[layer]])) {
			opt <- spatOptions()
			x@ptr <- x@ptr$makeCategorical(layer-1, opt)
			x <- messages(x, "as.factor<-")
		}
		if (is.data.frame(value)) {
			stopifnot(NCOL(value) == 2)
			x@ptr$setCategories(layer-1, value[,1], value[,2])
		} else if (is.vector(value)){
			x@ptr$setCategories(layer-1, 0:(length(value)-1), as.character(value))
		}
		messages(x, "levels<-")
	}
)

