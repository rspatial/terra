
if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}
setMethod ("area" , "SpatRaster",
	function (x, ...) {
		error("area", "this method was removed. Use expanse or cellSize")
	}
)


#if (!isGeneric("setCats")) { setGeneric("setCats", function(x, ...) standardGeneric("setCats")) }

#setMethod ("setCats" , "SpatRaster",
#	function (x, ...) {
#		warn("setCats", "this function will be removed. Please can use 'levels<-' or 'set.cats' instead")
#		set.cats(x, ...)
#	}
#)


#setMethod("costDistance", signature(x="SpatRaster"),
#	function(x, target=0, scale=1, maxiter=50, filename="", ...) {
#		warn("costDistance", "'costDistance' was renamed to 'costDist'. 'costDistance' will be removed in a #future version")
#		costDist(x, target=target, scale=scale, maxiter=maxiter, filename=filename, ...)
#	}
#)

setMethod("gridDistance", signature(x="SpatRaster"),
	function(x, ...) {
		error("gridDistance", "'gridDistance' was renamed to 'gridDist'")
		#. 'gridDistance' will be removed in a future version")
		#gridDist(x, target=target, scale=scale, maxiter=maxiter, filename=filename, ...) 
	}
)


setMethod("focalCor", signature(x="SpatRaster"),
	function(x, ...) {
		error("focalCor", "'focalCor' was renamed to 'focalPairs'")
		# focalPairs(x, ...)
	}
)
