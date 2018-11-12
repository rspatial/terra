# This file is part of the "terra" package.
# Copyright (c) 2018  Robert J. Hijmans 
# License GPL v3


setMethod('reclassify', signature(x='SpatRaster', rcl='ANY'), 
function(x, rcl, include.lowest=FALSE, right=TRUE, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...) {
	
	if ( is.null(dim(rcl)) ) { 
		stopifnot((length(rcl) %% 3 == 0))
		rcl <- matrix(rcl, ncol=3, byrow=TRUE) 
	} else if (is.data.frame(rcl)) {
		rcl <- as.matrix(rcl)
	}

	right <- ifelse(is.na(right), 2, ifelse(right, 0, 1))
	include.lowest <- as.logical(include.lowest[1])
	filename <- as.character(filename[1])

    x@ptr <- x@ptr$rcppReclassify(rcl, right, include.lowest, filename, format[1], datatype[1], overwrite[1])
	show_messages(x, "reclassify")	
}
)


		
