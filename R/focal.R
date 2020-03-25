# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, na.rm=TRUE, fillvalue=NA, fun="sum", filename="", overwrite=FALSE, wopt=list(), ...)  {

	cpp <- FALSE
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		w <- dim(w)
		na.rm <- FALSE
		fun <- "sum"
		cpp <- TRUE
	} else {
		w <- rep_len(w, 2)
		m <- 0.5[0]
		fun <- .makeTextFun(match.fun(fun))
		if (class(fun) == "character") { 
			test <- match(fun, c("mean", "min", "max", "sum", "median", "modal"))
			cpp <- TRUE
		}
	}
	
	if (cpp) {		
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$focal(w, m, fillvalue, na.rm, fun, opt)
		show_messages(x, "focal")
		return(x)
	
	} else {
		out <- rast(x)
		b <- writeStart(out, filename, overwrite, wopt)
		for (i in 1:b$n) {
			v <- matrix(x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i]), ncol=prod(w), byrow=TRUE)
			v <- apply(v, 1, fun, na.rm=na.rm)
			writeValues(out, v, c(b$row[i], b$nrows[i]))
		}
		out <- writeStop(out)		
		return(out)
	}
}
)
