# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, fun="sum", na.rm=TRUE, na.only=FALSE, fillvalue=NA, filename="", method=3, ...)  {

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

	if (is.matrix(w)) {
		m <- as.vector(t(w))
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		m <- rep(1, prod(w))
	}
	cpp <- FALSE
	txtfun <- .makeTextFun(match.fun(fun))
	if (is.character(txtfun)) { 
		if (txtfun %in% c("mean", "sum")) {
			cpp <- TRUE
		}
	}

	if (cpp) {
		opt <- spatOptions(filename, ...)
		if (method==1) {
			x@ptr <- x@ptr$focal1(w, m, fillvalue, na.rm[1], na.only[1], txtfun, opt)		
		} else if (method == 2) {
			x@ptr <- x@ptr$focal2(w, m, fillvalue, na.rm[1], na.only[1], txtfun, opt)
		} else {
			x@ptr <- x@ptr$focal3(w, m, fillvalue, na.rm[1], na.only[1], txtfun, opt)
		}
		messages(x, "focal")
		return(x)

	} else {
		out <- rast(x)
		readStart(x)
		on.exit(readStop(x))
		b <- writeStart(out, filename, ...)
		dow <- !all(m == 1)
		for (i in 1:b$n) {
			v <- matrix(x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i]), ncol=prod(w), byrow=TRUE)
			if (dow) {
				v <- v * w
			}
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
