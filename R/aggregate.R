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


setMethod('aggregate', signature(x='SpatRaster'), 
function(x, fact=2, fun='mean', na.rm=TRUE, filename="", overwrite=FALSE, wopt=list(), ...)  {

	#expand=TRUE, 
	fun <- .makeTextFun(match.fun(fun))
	toc <- FALSE
	if (class(fun) == "character") { 
		if (fun %in% c('sum', 'mean', 'min', 'max')) {
			toc <- TRUE
		}
	}

	opt <- .runOptions(filename[1], overwrite[1], wopt)	
	if (toc) {	
		#	fun='mean', expand=TRUE, na.rm=TRUE, filename=""
		x@ptr <- x@ptr$aggregate(fact, fun, na.rm, opt)
		return (show_messages(x, "aggregate"))
	} else {
		out <- rast(x)
		out@ptr <- out@ptr$aggregate(fact, "sum", na.rm, opt)
		out <- show_messages(out, "aggregate")
		dims <- x@ptr$get_aggregate_dims(fact)
		b <- x@ptr$getBlockSize(4)		
		nc <- ncol(x)	
		readStart(x)
		ignore <- writeStart(out, filename[1], overwrite[1], opt)
		for (i in 1:b$n) {
			v <- x@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
			v <- x@ptr$get_aggregates(v, b$nrows[i], dims)
			v <- sapply(v, fun, na.rm=na.rm)
			writeValues(out, v, b$row[i])
		}
		writeStop(out)
		readStop(x)
		return(out)
	}
}
)

