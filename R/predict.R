# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3

	

setMethod('predict', signature(object='SpatRaster'), 
	function(object, model, fun=predict, ..., filename="", overwrite=FALSE, wopt=list()) {
		out <- rast(object)
		readStart(object)
		nc <- ncol(object)
		b <- writeStart(out, filename, overwrite, wopt)
		for (i in 1:b$n) {
			d <- object@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
			d <- data.frame(matrix(d, ncol = ncol(object)))
			names(d) <- names(object)
			r <- predict(model, d, ...)
			writeValues(out, r, b$row[i])
		}
		writeStop(out)
		readStop(object)
		return(out)
	}
)


