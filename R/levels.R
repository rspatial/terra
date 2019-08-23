
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


setRat <- function(x, rat) {
	stopifnot(nlyr(x) == 1)
	#stopifnot(layer > 0 & layer <= nlyr(x))
	rat <- .makeSpatDF(rat)
	x@ptr$setAttributes(0, rat)
}	

	

setMethod("as.factor", signature(x="SpatRaster"), 
	function(x) {
		stopifnot(nlyr(x) == 1)
		stopifnot(hasValues(x))
		x@ptr <- x@ptr$copy()
		x@ptr$createCategories(0)
		show_messages(x)
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



setMethod("levels<-", signature(x="SpatRaster"), 
	function(x, value) {
		stopifnot(nlyr(x) == 1)
		#stopifnot(is.factor(x))
		#stopifnot(hasValues(x))
		if (is.data.frame(value)) {
			stopifnot(NCOL(value) == 2)
			x@ptr$setCategories(0, value[,1], value[,2])
		} else if (is.vector(value)){
			x@ptr$setCategories(0, 1:length(value), as.character(value))		
		}
		x <- show_messages(x)
	}
)

