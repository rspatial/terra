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
function(x, y, ..., fun, na.rm=TRUE, filename="", format="", datatype="FLT4S", overwrite=FALSE)  {

	out <- rast(x)
	readStart(x)
	b <- writeStart(out, filename, format, datatype, overwrite)
	for (i in 1:b$n) {
		v <- x@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
		vv <- as.list(data.frame(v))
		names(vv) <- NULL
		r <- do.call(fun, vv, na.rm=na.rm)
		# if i==1, check size of output and ajust layers
		writeValues(out, r, b$row[i])
	}
	writeStop(out)
	readStop(x)
	return(out)
}
)


setMethod("reduce", signature(x="SpatRaster"), 
function(x, fun, na.rm=TRUE, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...)  {
	
	txtfun <- .makeTextFun(match.fun(fun))
	if (class(txtfun) == 'character') { 
		if (txtfun %in% c("max", "min", "range", "prod", "sum", "any", "all"))
		x@ptr <- x@ptr$summary(txtfun, na.rm, filename, format, datatype, overwrite)	
		return(x)
	} 
		
	out <- rast(x)
	readStart(x)
	b <- writeStart(out, filename, format, datatype, overwrite)
	for (i in 1:b$n) {
		v <- x@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
		r <- apply(v, 1, fun, na.rm=na.rm)
		# if i==1, check size of output and ajust layers
		writeValues(out, r, b$row[i])
	}
	writeStop(out)
	readStop(x)
	return(out)
}
)
