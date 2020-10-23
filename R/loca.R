
setMethod("loca", signature(x="SpatRaster"), 
function(x, fun, ..., nodes=1, filename="", overwrite=FALSE, wopt=list())  {

	out <- rast(x)
	nc <- ncol(x)
	readStart(x)
	on.exit(readStop(x))

# test the shape of the output by testing with one row
	v <- readValues(x, round(0.5*nrow(x)), 1, 1, nc, mat=TRUE)
	r <- try(fun(as.vector(v), ...))
	if (inherits(r, "try-error")) {
		stop("'fun' is not valid")
	}
	if (!is.vector(r)) {
		stop("'fun' does not return a vector")
	}
	if (!(is.numeric(r) | is.logical(r))) {
		stop("'fun' does not return a numeric vector")
	} 
	if (length(r) != length(v)) {
		stop("'fun' does not return the same number of values as the input")
	}

	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
		v <- fun(as.vector(v), ...)
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	writeStop(out)
}
)

