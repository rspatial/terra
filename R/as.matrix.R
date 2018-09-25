# Author: Robert J. Hijmans 
# Date : October 2018
# Version 1.0
# Licence GPL v3



# mode argument is ignored as mode=mode gave an error on R-devel
setMethod('as.vector', signature(x='SpatExtent'), 
function(x, mode='any') {
	x@ptr$vector
})

setMethod('as.character', signature(x='SpatExtent'), 
function(x, mode='any') {
	paste("'", x@ptr$vector, collapse="', '", "'")
})


setMethod('as.vector', signature(x='SpatRaster'), 
function(x, mode='any') {
	values(x, FALSE)
})



setMethod('as.matrix', signature(x='SpatRaster'), 
function(x, mode='any') {
	matrix( values(x, FALSE), nrow=nrow(x), byrow=TRUE)
})


setMethod('as.array', signature(x='SpatRaster'), 
function(x, mode='any') {
	dm <- dim(x)
	x <- values(x, TRUE)
	ar <- array(NA, dm)
	for (i in 1:dm[3]) {
		ar[,,i] <- matrix(x[,i], nrow=dm[1], byrow=TRUE)
	}
	ar		
})


