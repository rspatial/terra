# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# Licence GPL v3


setMethod("c", signature(x="SpatRaster"), 
	function(x, ...) {
		for (i in list(...)) {
			if (class(i) == "SpatRaster") {
				x@ptr <- x@ptr$combineSources(i@ptr)
			}
		}
		x@ptr$setNames(x@ptr$names)
		show_messages(x, "c")		
	}
)


setMethod("clamp", signature(x="SpatRaster"), 
	function(x, lower=-Inf, upper=Inf, values=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename[1], overwrite[1],wopt)
		x@ptr <- x@ptr$clamp(lower, upper, values[1], opt)
		show_messages(x, "clamp")
	}
)

setMethod("crop", signature(x="SpatRaster", y="ANY"), 
	function(x, y, snap="near", filename="", overwrite=FALSE, wopt=list(), ...) {
		if (!inherits(y, "SpatExtent")) {
			y <- try(ext(y))
			if (class(y) == "try-error") { stop("cannot get an extent from y") }
		}
		opt <- .runOptions(filename[1], overwrite[1],wopt)
		x@ptr <- x@ptr$crop(y@ptr, snap[1], opt)
		show_messages(x, "crop")		
	}
)


setMethod("gridDistance", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename[1], overwrite[1],wopt)
		x@ptr <- x@ptr$gridDistance(opt)
		show_messages(x, "gridDistance")
	}
)


setMethod("mask", signature(x="SpatRaster", mask="SpatRaster"), 
	function(x, mask, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename[1], overwrite[1],wopt)
		x@ptr <- x@ptr$mask(mask@ptr, opt)
		show_messages(x, "mask")		
	}
)

setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, background=NA, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename[1], overwrite[1],wopt)
		y@ptr <- y@ptr$rasterizePolygons(x@ptr, background[1], opt)
		show_messages(y, "rasterize")
	}
)


setMethod('reclassify', signature(x='SpatRaster', rcl='ANY'), 
function(x, rcl, include.lowest=FALSE, right=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
	
	
	if ( is.null(dim(rcl)) ) { 
		stopifnot((length(rcl) %% 3 == 0))
		rcl <- matrix(rcl, ncol=3, byrow=TRUE) 
	} else if (is.data.frame(rcl)) {
		rcl <- as.matrix(rcl)
	}

	right <- ifelse(is.na(right), 2, ifelse(right, 0, 1))
	include.lowest <- as.logical(include.lowest[1])

	opt <- .runOptions(filename[1], overwrite[1],wopt)
    x@ptr <- x@ptr$rcppReclassify(rcl, right, include.lowest, opt)
	show_messages(x, "reclassify")	
}
)


setMethod("sampleRegular", signature(x="SpatRaster", size="numeric"), 
	function(x, size, ...) { 
		x@ptr <- x@ptr$sampleRegular(size)
		show_messages(x, "sampleRegular")		
	}
)


setMethod("trim", signature(x="SpatRaster"), 
	function(x, padding=0, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename[1], overwrite[1], wopt)
		x@ptr <- x@ptr$trim(padding[1], opt)
		show_messages(x, "rasterize")
	}
)


