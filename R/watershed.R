# Author: Emanuele Cordano
# Date : October 2023
# Version 1.0
# License GPL v3

setMethod("watershed", signature(x="SpatRaster"), 
    function(x, pourpoint, filename="", ...) { 
        opt <- spatOptions(filename, ...)		
		cell <- cellFromXY(x, pourpoint)
		if (is.na(cell)) error("watershed", "pourpoint not on raster")
        x@pntr <- x@pntr$watershed2(as.integer(cell-1), opt)
        messages(x, "watershed") ## EC 20210318
    }
)

setMethod("pitfinder", signature(x="SpatRaster"), 
    function(x, filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@pntr <- x@pntr$pitfinder2(opt)
        messages(x, "pitfinder") ## EC 20210318
    }
)

setMethod("NIDP", signature(x="SpatRaster"), 
    function(x, filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@pntr <- x@pntr$NIDP2(opt)
        messages(x, "NIDP") ## EC 20231031
    }
)

setMethod("flowAccumulation", signature(x="SpatRaster"), 
    function(x, weight=NULL, filename="", ...) { 
		opt <- spatOptions(filename, ...)
		if (is.null(weight)) {      
			x@pntr <- x@pntr$flowAccu2(opt)
	    } else {
			x@pntr <- x@pntr$flowAccu2_weight(weight@pntr, opt)
		} 
		messages(x, "flowAccumulation") 
    }      
)


