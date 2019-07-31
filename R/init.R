# Author: Robert J. Hijmans
# Date :  May 2019
# Version 1.0
# License GPL v3

	
setMethod("init", signature(x="SpatRaster"), 
	function(x, fun, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		if (is.character(fun)) {
			fun <- fun[1]
			if (fun %in% c("x", "y", "row", "col", "cell", "chess")) {
				x@ptr <- x@ptr$initf(fun, TRUE, opt)
				show_messages(x, "init")
			} else {
				stop("unknown function")
			}
		} else if (is.numeric(fun)) {
			x@ptr <- x@ptr$initv(fun[1], opt)
			show_messages(x, "init")
		} else {
			out <- rast(x)
			nc <- ncol(out)
			b <- writeStart(out, filename, overwrite, wopt)
			for (i in 1:b$n) {
				n <- b$nrows[i] * nc;
				r <- fun(n)
				writeValues(out, r, b$row[i])
			}
			out <- writeStop(out)
			return(out)		
		}
	}
)
	
