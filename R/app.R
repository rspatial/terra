
setMethod("app", signature(x="SpatRaster"), 
function(x, fun, ..., filename="", overwrite=FALSE, wopt=list())  {


	txtfun <- .makeTextFun(match.fun(fun))
	if (class(txtfun) == "character") { 
		if (txtfun %in% c("max", "min", "mean", "range", "prod", "sum", "any", "all")) {
			opt <- .runOptions(filename, overwrite, wopt)
			na.rm <- isTRUE(list(...)$na.rm)
			x@ptr <- x@ptr$summary(txtfun, na.rm, opt)	
			return(show_messages(x, "app"))
		}		
	}

	out <- rast(x)
	nlyr(out) <- 1
	nc <- ncol(x)
	readStart(x)

# figure out the shape of the output by testing with one row
	v <- readValues(x, round(0.5*nrow(x)), 1, 1, nc, mat=TRUE)
	#narg <- sum(sapply(f, as.character) == "", na.rm=TRUE)
	#if (narg > 1) {
	#	vv <- as.list(as.data.frame(v))
	#	r <- do.call(fun, vv, ...)	
	#} else {
	r <- apply(v, 1, fun, ...)
	
	#}
	if (is.list(r)) {
		stop("the function returns a list (should be numeric or matrix")
	}
	trans <- FALSE			
	if (NCOL(r) > 1) {
		#? if ((ncol(r) %% nc) == 0) {
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
		writeValues(out, r, b$row[i], b$nrows[i])
	}
	readStop(x)
	writeStop(out)
}
)


