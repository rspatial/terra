# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, na.rm=TRUE, na.only=FALSE, fillvalue=NA, fun="sum", filename="", overwrite=FALSE, wopt=list(), ...)  {

	if (nlyr(x) > 1) {
		warn("focal", "only the first layer of x is used")
		x <- x[[1]]
	}
	if (na.only && (!is.matrix(w))) {
		if (!na.rm) {
			warning("focal", "na.rm set to TRUE because na.only is TRUE")
			na.rm <- TRUE
		}
	}

	cpp <- FALSE
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		w <- dim(w)
		#na.rm <- FALSE
		fun <- .makeTextFun(match.fun(fun))
		if (!is.character(fun))	fun <- "bad"
		if (fun != "sum") {
			warning("focal", "if 'w' is a matrix, 'fun' must be 'sum'")
		}
		fun <- "sum"
		cpp <- TRUE
	} else {
		w <- rep_len(w, 2)
		m <- 0.5[0]
		fun <- .makeTextFun(match.fun(fun))
		if (class(fun) == "character") { 
			test <- match(fun, c("mean", "min", "max", "sum", "median", "modal", "sd", "sdpop"))
			cpp <- TRUE
		}
	}

	if (cpp) {
		opt <- spatOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$focal(w, m, fillvalue, na.rm[1], na.only[1], fun, opt)
		messages(x, "focal")
		return(x)

	} else {
		out <- rast(x)
		readStart(x)
		on.exit(readStop(x))
		b <- writeStart(out, filename, overwrite, wopt)
		for (i in 1:b$n) {
			v <- matrix(x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i]), ncol=prod(w), byrow=TRUE)
			v <- apply(v, 1, fun, na.rm=na.rm)
			if (na.only) {
				vv <- readValues(x, b$row[i], b$nrows[i])
				j <- !is.na(vv)
				v[j] <- vv[j]
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		}
		out <- writeStop(out)
		return(out)
	}
}
)
