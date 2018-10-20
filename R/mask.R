# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# Licence GPL v3


setMethod('mask', signature(x='SpatRaster', mask='SpatRaster'), 
function(x, mask, filename="", overwrite=FALSE, ...) { 
    x@ptr <- x@ptr$mask(mask@ptr, filename[1], overwrite[1])
	.messages(x, "mask")		
	x	
}
)

