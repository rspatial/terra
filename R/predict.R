# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3

	
setMethod("predict", signature(object="SpatRaster"), 
	function(object, model, fun=predict, nlyr, ..., filename="", overwrite=FALSE, wopt=list()) {
		
		if (missing(nlyr)) {
			nlyr <- try(NCOL(fun(model, object[1:min(10,ncell(object))])), silent=TRUE)
			if (class(nlyr) == "try-error") {
				nlyr=1
			} else {
				nlyr = max(1, nlyr)
			}
		}
		
		out <- rast(object, nlyr=nlyr)
		readStart(object)
		nc <- ncol(object)
		b <- writeStart(out, filename, overwrite, wopt)
		on.exit(writeStop(out))		
		for (i in 1:b$n) {
			d <- object@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
			d <- data.frame(matrix(d, ncol = nlyr(object)))
			names(d) <- names(object)
			r <- fun(model, d, ...)
			writeValues(out, r, b$row[i])
		}
		readStop(object)
		return(out)
	}
)


