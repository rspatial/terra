# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3

setMethod("as.list", signature(x="SpatRaster"), 
function(x, ...)  {
	for (i in 1:nlyr(x)) {
		out[[i]] <- x[[i]]
	}
	out
}
)



setMethod("overlay", signature(x="SpatRaster", y="SpatRaster"), 
function(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())  {
	
	stopifnot(!missing(fun))
	
	compareGeom(x, y, lyrs=FALSE, crs=TRUE, warncrs=TRUE, ext=TRUE, rowcol=TRUE, res=FALSE)
	
	nl <- max(nlyr(x), nlyr(y))
	out <- rast(x,nlyr=nl)

	readStart(x)
	readStart(y)
	b <- writeStart(out, filename, overwrite, wopt)
	on.exit(writeStop(out))
	nc <- ncol(x)
	nl <- nlyr(x)
#	fnames <- names(formals(fun))
#	if (length(fnames) != nl) {	dnames <- NULL 	} else { dnames <- list(list(), fnames)	}

	for (i in 1:b$n) {
		vx <- x@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
		vy <- y@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
		r <- fun(vx, vy)
		writeValues(out, r, b$row[i])
	}
	readStop(x)
	return(out)
}
)

