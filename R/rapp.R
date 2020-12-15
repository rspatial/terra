
setMethod("rapp", signature(x="SpatRaster"), 
function(x, index, fun, ..., filename="", overwrite=FALSE, wopt=list()) {

	stopifnot(inherits(index, "SpatRaster"))
	txtfun <- .makeTextFun(match.fun(fun))
	if (inherits(txtfun, "character")) { 
		if (txtfun %in% c("max", "min", "mean", "prod", "sum", "any", "all")) {
			opt <- spatOptions(filename, overwrite, wopt)
			na.rm <- isTRUE(list(...)$na.rm)
			x@ptr <- x@ptr$rapply(index@ptr, txtfun, na.rm, opt)	
			return(messages(x, "rapp"))
		}		
	} 

	stopifnot(hasValues(x))
	stopifnot(hasValues(index))
	stopifnot(nlyr(index) == 2)
	compareGeom(x, index, lyrs=FALSE, crs=FALSE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE) 

	out <- rast(x)
	nlyr(out) <- 1
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- x@ptr$rappvals(index@ptr, b$row[i]-1, b$nrows[i])
		v <- sapply(v, fun, ...)
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)


