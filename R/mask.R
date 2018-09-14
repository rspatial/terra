# Author: Robert J. Hijmans
# Date : November 2009
# Version 1.0
# Licence GPL v3


setMethod('mask', signature(x='SpatRaster', mask='SpatRaster'), 
function(x, mask, filename="", ...) { 
    overwrite <- .overwrite(...)
    r <- methods::new("SpatRaster")
    r@ptr <- x@ptr$mask(mask@ptr, filename[1], overwrite)
    return(r)
}
)

