# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# Licence GPL v3


setMethod('c', signature(x='SpatRaster'), 
	function(x, ...) {
		for (i in list(...)) {
			if (class(i) == 'SpatRaster') {
				x@ptr <- x@ptr$combineSources(i@ptr)
			}
		}
		x@ptr$setNames(x@ptr$names)
		show_messages(x, "c")		
	}
)


setMethod("crop", signature(x="SpatRaster", y="ANY"), 
function(x, y, snap="near", filename="", format="", datatype="FLT4S", overwrite=FALSE, ...) {
	if (!inherits(y, "SpatExtent")) {
		y <- try(ext(y))
		if (class(y) == "try-error") { stop("cannot get an extent from y") }
	}
	x@ptr <- x@ptr$crop(y@ptr, snap[1], filename[1], format[1], datatype[1], overwrite[1])
	show_messages(x, "crop")		
}
)


setMethod("mask", signature(x="SpatRaster", mask="SpatRaster"), 
function(x, mask, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...) { 
    x@ptr <- x@ptr$mask(mask@ptr, filename[1], format[1], datatype[1], overwrite[1])
	show_messages(x, "mask")		
}
)

setMethod("rasterize", signature(x="SpatLayer", y="SpatRaster"), 
function(x, y, background=NA, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...) { 
	y@ptr <- y@ptr$rasterizePolygons(x@ptr, background[1], filename[1], format[1], datatype[1], overwrite[1])
	show_messages(y, "rasterize")
}
)

setMethod("sampleRegular", signature(x="SpatRaster", size="numeric"), 
function(x, size, ...) { 
    x@ptr <- x@ptr$sampleRegular(size)
	show_messages(x, "sampleRegular")		
}
)


setMethod('trim', signature(x='SpatRaster'), 
function(x, padding=0, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...) {
	x@ptr <- x@ptr$trim(padding[1], filename[1], format[1], datatype[1], overwrite[1])
	show_messages(x, "rasterize")
}
)


