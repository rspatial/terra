# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# Licence GPL v3


.makeTextFun <- function(fun) {
	if (class(fun) != 'character') {
		if (is.primitive(fun)) {
			test <- try(deparse(fun)[[1]], silent=TRUE)
			if (test == '.Primitive(\"sum\")') { fun <- 'sum' 
			} else if (test == '.Primitive(\"min\")') { fun <- 'min' 
			} else if (test == '.Primitive(\"max\")') { fun <- 'max' 
			}
		} else {
			test1 <- isTRUE(try( deparse(fun)[2] == 'UseMethod(\"mean\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == 'mean', silent=TRUE))
			if (test1 | test2) { 
				fun <- 'mean' 
			}
		} 
	}
	return(fun)
}

.overwrite <- function(...) {
	overwrite <- list(...)$overwrite
	if (is.null(overwrite)) { 
		overwrite <- FALSE 
	}
	overwrite
}

setMethod('aggregate', signature(x='SpatRaster'), 
function(x, fact=2, fun='mean', na.rm=TRUE, filename="", ...)  {

	#expand=TRUE, 
	overwrite <- .overwrite(...)
	
	fact <- round(fact)
	lf <- length(fact)
	if (lf > 3) {
		stop('fact should have length 1, 2, or 3')
	}
	if (any(fact < 1)) {
		stop('fact should be > 0')
	}
	if (! any(fact > 1)) {
		warning('all fact(s) were 1, nothing to aggregate')
		return(x)
	}
	dims <- x@ptr$get_aggregate_dims(fact)
	fun <- .makeTextFun(match.fun(fun))
	if (class(fun) == 'character') { 
		op <- as.integer(match(fun, c('sum', 'mean', 'min', 'max')) - 1)
	} else {
		op <- NA
	}

	if (!is.na(op)) {	
		r <- methods::new('SpatRaster')
		#	fun='mean', expand=TRUE, na.rm=TRUE, filename=""
		x@ptr <- x@ptr$aggregate(dims, fun, na.rm, filename, overwrite)
		return (show_messages(x, "aggregate"))
	} else {
		e <- as.vector(ext(x))
		rs <- res(x)
		e[2] <- e[1] + dims[4] * rs[2];
		e[3] <- e[4] - dims[5] * rs[1];
		out <- rast(nrow=dims[4], ncol=dims[4], nlyr=dims[6],  crs=crs(x), ext=e)
		nc <- ncol(x)
		
		readStart(x)
		b <- writeStart(out, filename)
		for (i in 1:b$n) {
			#v <- x@ptr$get_aggregates(dims, b$row[i], b$nrows[i], 1, nc)
			v <- x@ptr$get_aggregates(dims)
			v <- do.call(rbind, v)
			v <- apply(v, 1, fun, na.rm=na.rm)
			writeValues(out, v, b$row[i])
		}
		writeStop(out)
		readStop(x)
		return(out)
	}
}
)

