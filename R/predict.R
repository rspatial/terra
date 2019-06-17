# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3

	
setMethod("predict", signature(object="SpatRaster"), 
	function(object, model, fun=predict, ..., filename="", overwrite=FALSE, wopt=list()) {
		
		#if (missing(nlyr)) {
		#	nlyr <- try(NCOL(fun(model, object[1:min(10,ncell(object))])), silent=TRUE)
		#	if (class(nlyr) == "try-error") {
		#		nlyr = 1
		#	} else {
		#		nlyr = max(1, nlyr)
		#	}
		#}
		
		nl <- 1
		nc <- ncol(object)
		tomat <- FALSE
		v <- readValues(object, round(0.5*nrow(object)), 1, 1, nc, TRUE, TRUE)
		v <- fun(model, v, ...)
		if (NCOL(v) > 1) {
			nl <- ncol(v)
			if (inherits(v, "data.frame")) {
				tomat <- TRUE
			}
		}
		
		out <- rast(object, nlyr=nl)
		readStart(object)
		b <- writeStart(out, filename, overwrite, wopt)
		for (i in 1:b$n) {
			d <- readValues(object, b$row[i], b$nrows[i], 1, nc, TRUE, TRUE)
			r <- fun(model, d, ...)
			if (tomat) {
				r <- as.matrix(r)
			}
			writeValues(out, r, b$row[i])
		}
		readStop(object)
		writeStop(out)
		return(out)
	}
)


