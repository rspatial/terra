
setMethod("reduce", signature(x="SpatRaster", fun="function"), 
function(x, fun, ..., filename="", overwrite=FALSE, wopt=list())  {

	opt <- terra:::.runOptions(filename, overwrite, wopt)

	txtfun <- terra:::.makeTextFun(match.fun(fun))
	if (class(txtfun) == "character") { 
		if (txtfun %in% c("max", "min", "mean", "range", "prod", "sum", "any", "all")) {
			narm <- ifelse(isTRUE(list(...)$na.rm), TRUE, FALSE)
			x@ptr <- x@ptr$summary(txtfun, narm, opt)	
			return(show_messages(x))
		}		
	}

	out <- rast(x)
	nlyr(out) <- 1
	readStart(x)
	nc <- ncol(x)

# figure out the shape of the output by testing with one row
	v <- readValues(x, round(0.5*nrow(x)), 1, 1, nc, TRUE)
	r <- apply(v, 1, fun, ...)
	trans <- FALSE			
	if (NCOL(r) > 1) {
		if (ncol(r) == nc) {
			nlyr(out) <- nrow(r)
			trans <- TRUE
		} else if (nrow(r) == nc) {
			nlyr(out) <- ncol(r)
		} else {
			stop("cannot handle this function")
		}
	}

	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
		r <- apply(v, 1, fun, ...)
		if (trans) {
			r <- t(r)
			#r <- as.vector(r)
		}
		writeValues(out, r, b$row[i])
	}
	writeStop(out)
	readStop(x)
	return(out)
}
)


