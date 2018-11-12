# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# Licence GPL v3

setMethod("as.list", signature(x="SpatRaster"), 
function(x, ...)  {
	for (i in 1:nlyr(x)) {
		out[[i]] <- x[[i]]
	}
	out
}
)


setMethod("overlay", signature(x="SpatRaster", y="SpatRaster"), 
function(x, y, fun, ..., filename="", format="", datatype="FLT4S", overwrite=FALSE)  {
	
	stopifnot(!missing(fun))
	x@ptr$compare_geom(y@ptr, FALSE, TRUE, TRUE)
	x <- show_messages(x)
	
	nl <- max(nlyr(x), nlyr(y))
	out <- rast(x,nlyr=nl)

	readStart(x)
	readStart(y)
	b <- writeStart(out, filename, format, datatype, overwrite)
#	b <- writeStart(out, "", "", "", FALSE)
	nc <- ncol(x)
	nl <- nlyr(x)

#	fnames <- names(formals(fun))
#	if (length(fnames) != nl) {	dnames <- NULL 	} else { dnames <- list(list(), fnames)	}

	for (i in 1:b$n) {
		vx <- x@ptr$readValues(b$row[i]-1, b$nrows[i], 0, nc)
		vy <- y@ptr$readValues(b$row[i]-1, b$nrows[i], 0, nc)
		r <- fun(vx, vy)
		writeValues(out, r, b$row[i])
	}
	writeStop(out)
	readStop(x)
	return(out)
}
)

