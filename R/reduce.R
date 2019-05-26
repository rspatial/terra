
setMethod("reduce", signature(x="SpatRaster", fun="function"), 
function(x, fun, ..., filename="", overwrite=FALSE, wopt=list())  {

	opt <- .runOptions(filename, overwrite,wopt)

	txtfun <- terra:::.makeTextFun(match.fun(fun))
	if (class(txtfun) == 'character') { 
		if (txtfun %in% c("max", "min", "range", "prod", "sum", "any", "all"))
		na.rm = ifelse(isTRUE(list(...)$na.rm), TRUE, FALSE)
		x@ptr <- x@ptr$summary(txtfun, na.rm, opt)	
		return(x)
	} 
	

	out <- rast(x)
	readStart(x)
	nc <- ncol(x)
	b <- writeStart(out, opt)
	for (i in 1:b$n) {
		v <- x@ptr$readValues(b$row[i], b$nrows[i], 0, nc)
		r <- apply(v, 1, fun, ...)
		# if i==1, check size of output and adjust layers
		writeValues(out, r, b$row[i])
	}
	writeStop(out)
	readStop(x)
	return(out)
}
)


