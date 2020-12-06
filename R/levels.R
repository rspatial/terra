
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


setRat <- function(x, layer, rat) {
	#stopifnot(nlyr(x) == 1)
	stopifnot(layer > 0 & layer <= nlyr(x))
	rat <- .makeSpatDF(rat)
	x@ptr$setAttributes(layer-1, rat)
}	

	

setMethod("as.factor", signature(x="SpatRaster"), 
	function(x) {
		stopifnot(nlyr(x) == 1)
		stopifnot(hasValues(x))
		x@ptr <- x@ptr$deepcopy()
		x@ptr$createCategories(0)
		show_messages(x, "as.factor")
	}
)

setMethod("is.factor", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasCategories()
	}
)


setMethod("levels", signature(x="SpatRaster"), 
	function(x) {
		#if (any(x@ptr$hasCategories())) {
		x@ptr$getCategories()
		#} else {
		#	vector(mode="list", nlyr(x))
		#}
	}
)


setMethod("levels<-", signature(x="SpatRaster"), 
	function(x, value) {
		stopifnot(nlyr(x) == 1)
		if (is.null(value) | is.na(value[1])) {
			x@ptr$removeCategories(0)
		}
		#stopifnot(is.factor(x))
		#stopifnot(hasValues(x))
		if (is.data.frame(value)) {
			stopifnot(NCOL(value) == 2)
			x@ptr$setCategories(0, value[,1], value[,2])
		} else if (is.vector(value)){
			x@ptr$setCategories(0, as.character(value), 0:(length(value)-1))		
		}
		x <- show_messages(x, "levels<-")
	}
)

