
setMethod("as.factor", signature(x="SpatRaster"), 
	function(x) {
		stopifnot(nlyr(x) == 1)
		stopifnot(hasValues(x))
		x@ptr <- x@ptr$copy()
		x@ptr$createCategories(0)
		show_messages(x)
	}
)

setMethod('is.factor', signature(x='SpatRaster'), 
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

