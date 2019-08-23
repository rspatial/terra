


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


setRat <- function(x, layer=0, rat) {
	stopifnot(layer > 0 & layer <= nlyr(x))
	rat <- .makeSpatDF(rat)
	x@ptr$setAttributes(layer, rat)
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
		stopifnot(is.factor(x))
		stopifnot(hasValues(x))
		x@ptr$setCategories(0, value)
		x <- show_messages(x)
	}
)

