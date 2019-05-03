# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, na.rm=TRUE, fillvalue=NA, fun="sum", filename="", overwrite=FALSE, wopt=list(), ...)  {

	if (is.matrix(w)) {
		stopifnot(ncol(w) == nrow(w))	
		w <- as.vector(t(w))  # row-wise
		cpp <- TRUE
		fun <- "sum"
	} else {
		w <- rep_len(1, w^2)
		fun <- .makeTextFun(match.fun(fun))
		if (class(fun) == "character") { 
			test <- match(fun, c("mean", "min", "max", "sum"))
			cpp <- TRUE
		} else {
			cpp <- FALSE
		}
	}
	

	if (cpp) {
		
		opt <- .runOptions(filename[1], overwrite[1], wopt)
		x@ptr <- x@ptr$focal(w, fillvalue, na.rm, fun, opt);
		show_messages(x, "focal")
		return(x)
	
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
)
