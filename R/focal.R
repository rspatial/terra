# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod('focal', signature(x='SpatRaster'), 
function(x, w=3, na.rm=TRUE, fillvalue=NA, fun='sum', filename="", overwrite=FALSE, wopt=list(), ...)  {

	if (is.matrix(w)) {
		stopifnot(ncol(w) == nrow(w))	
		w <- as.vector(t(w))  # row-wise
		domat <- TRUE
	} else {
		w <- rep_len(w, 3)
		fun <- .makeTextFun(match.fun(fun))
		if (class(fun) == 'character') { 
			op <- as.integer(match(fun, c('mean', 'min', 'max', 'sum')) - 1)
		} else {
			stop()
		}
		domat <- FALSE		
	}
	
	if (domat) {
		r <- methods::new('SpatRaster')
		ptr <- try(x@ptr$focal(w, fillvalue, na.rm, fun=3, filename, overwrite));
		if (class(ptr) == 'try-error') { stop("focal error") } else { r@ptr <- ptr }
		return(r)
	
	} else {
		if (!is.na(op)) {	
			r <- methods::new('SpatRaster')
			#	fun='mean', expand=TRUE, na.rm=TRUE, filename=""
			ptr <- try(x@ptr$focal(w, fillvalue, na.rm, fun, filename, overwrite));
			if (class(ptr) == 'try-error') { stop("focal error") } else { r@ptr <- ptr }
			return(r)
		} else {
			out <- rast(x)
			readStart(x)
			b <- writeStart(out, filename, overwrite, wopt)
			for (i in 1:b$n) {
				v <- matrix(x@ptr$focalValues(b$row[i], b$nrows[i], w), ncol=prod(w), byrow=TRUE)
				v <- apply(v, 1, fun, na.rm=na.rm)
				writeValues(out, v, b$row[i])
			}
			writeStop(out)
			readStop(x)
			return(out)
		}
	}
}	
)
