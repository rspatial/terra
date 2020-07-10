
if (!isGeneric("extend")) {setGeneric("extend", function(x, y, ...) standardGeneric("extend"))}	
setMethod("extend", signature(x="SpatRaster"), 
	function(x, y, ...) {
		stop("terra::extend has been removed. Use 'expand' instead")
	}
)

setMethod("expand", signature(x="SpatExtent"), 
function(x, y, ...) {
	if (length(y) == 1) {
		y <- rep(y, 4)
	} else if (length(y) == 2) {
		y <- rep(y, each=2)
	} else if (! length(y) == 4 ) {
		stop('argument "y" should be a vector of 1, 2, or 4 elements')	
	}
	e <- as.vector(x)
	e[1] <- e[1] - y[1]
	e[2] <- e[2] + y[2]
	e[3] <- e[3] - y[3]
	e[4] <- e[4] + y[4]
	ext(e)
}
)



setMethod('expand', signature(x='SpatRaster'), 
function(x, y, filename="", overwrite=FALSE, wopt=list(), ...) {

	if (!inherits(y, "SpatExtent")) {

		if (is.vector(y)) {
			if (length(y) <= 2) {
				adj <- abs(y) * rev(res(x))
				y <- as.vector(ext(x))
				y[1] <- y[1] - adj[1]
				y[2] <- y[2] + adj[1]
				y[3] <- y[3] - adj[2]
				y[4] <- y[4] + adj[2]
				y <- ext(y)
			}
		} else {
			test <- try ( y <- ext(y), silent=TRUE )
			if (class(test) == "try-error") {
				stop('Cannot get a SpatExtent object from argument y')
			}
		}
	}
	
	opt <- .runOptions(filename, overwrite, wopt)
	x@ptr <- x@ptr$expand(y@ptr, opt)
	show_messages(x, "expand")
}
)


