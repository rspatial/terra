# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# Licence GPL v3


setMethod('c', signature(x='SpatRaster'), 
	function(x, ...) {
		for (i in list(...)) {
			if (class(i) == 'SpatRaster') {
				x@ptr <- x@ptr$addSources(i@ptr)
			}
		}
		x@ptr$setNames(x@ptr$names)
		.messages(x, "c")		
	}
)


setMethod("crop", signature(x="SpatRaster", y="ANY"), 
function(x, y, snap="near", filename="", overwrite=FALSE, ...) {

	if (!inherits(y, "SpatExtent")) {
		y <- try(ext(y))
		if (class(y) == "try-error") { stop("cannot get an extent from y") }
	}
	x@ptr <- x@ptr$crop(y@ptr, filename, snap, overwrite)
	.messages(x, "crop")		
}
)

setMethod("extract", signature(x="SpatRaster", y="SpatLayer"), 
function(x, y, fun="", ...) { 
    r <- x@ptr$extract(y@ptr, fun)
	x <- .messages(x, "extract")		
	r
}
)


setMethod("mask", signature(x="SpatRaster", mask="SpatRaster"), 
function(x, mask, filename="", overwrite=FALSE, ...) { 
    x@ptr <- x@ptr$mask(mask@ptr, filename[1], overwrite[1])
	.messages(x, "mask")		
}
)

setMethod("rasterize", signature(x="SpatLayer", y="SpatRaster"), 
function(x, y, background=NA, filename="", ...) { 
	overwrite <- .overwrite(...)
	y@ptr <- y@ptr$rasterizePolygons(x@ptr, background, filename[1], overwrite)
	.messages(y, "rasterize")
}
)

setMethod('trim', signature(x='SpatRaster'), 
function(x, padding=0, filename='', ...) {
	overwrite <- .overwrite(...)
	x@ptr <- x@ptr$trim(padding, filename, overwrite)
	.messages(x, "rasterize")
}
)


