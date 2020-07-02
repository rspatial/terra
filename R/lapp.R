# Author: Robert J. Hijmans
# Date : May 2020
# Version 1.0
# License GPL v3

.lapp_test <- function(v, fun, usenames, ...) {
# figure out the shape of the output
	nms = ""
	nr <- nrow(v)
	if (!usenames) colnames(v) <- NULL
	vtst <- try(do.call(fun, c(v, list(...))), silent=FALSE)
	if (inherits(vtst, "try-error")) {
		nl <- -1
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {	
			nl <- -1
		}
	} else {
		nl <- -1
	}
	if (is.matrix(vtst)) {
		nms <- colnames(vtst)
	}
	list(nl=nl, names=nms)
}


setMethod("lapp", signature(x="SpatRaster"), 
function(x, fun, ..., usenames=FALSE, filename="", overwrite=FALSE, wopt=list())  {
	
	fun <- match.fun(fun)
	dots <- list(...)
	if (any(sapply(dots, function(i) inherits(i, "SpatRaster")))) {
		stop("Only 'x' can be a SpatRaster" )
		# otherwise .lapp_test may crash! 
	}
	
	if (usenames) {
		fnames <- names(formals(fun))
		x <- x[[names(x) %in% fnames]]
	}
	readStart(x)
	on.exit(readStop(x))
	ncx <- ncol(x)
	v <- readValues(x, round(0.5*nrow(x)), 1, 1, ncx, dataframe=TRUE)
	test <- .lapp_test(v, fun, usenames, ...)
	if (test$nl < 1) stop("lapp does not like 'fun'")
	out <- rast(x, nlyr=test$nl)
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, dataframe=TRUE)
		if (!usenames) colnames(v) <- NULL
		v <- do.call(fun, c(v, list(...)))
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)


.lapp_test_stack <- function(v, fun, recycle, ...) {
# figure out the shape of the output
	nms = ""
	nr <- nrow(v[[1]])
	if (recycle) {
		v <- lapply(v, as.vector)
	}
	vtst <- try(do.call(fun, c(v, list(...))), silent=FALSE)
	if (inherits(vtst, "try-error")) {
		nl <- -1
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {	
			nl <- -1
		}
	} else {
		nl <- -1
	}
	if (is.matrix(vtst)) {
		nms <- colnames(vtst)
	}
	list(nl=nl, names=nms)
}



setMethod("lapp", signature(x="SpatDataSet"), 
function(x, fun, ..., recycle=FALSE, filename="", overwrite=FALSE, wopt=list())  {
	
	fun <- match.fun(fun)
	dots <- list(...)
	if (any(sapply(dots, function(i) inherits(i, "SpatDataSet")))) {
		stop("Only 'x' can be a SpatDataSet" )
		# otherwise .lapp_test_stack fails
	}
	
	ncx <- ncol(x[1])
	nrx <- nrow(x[1])
	readStart(x)
	on.exit(readStop(x))
	
	v <- lapply(1:length(x), function(i) readValues(x[i], round(0.5*nrx), 1, 1, ncx, mat=TRUE))
	test <- .lapp_test_stack(v, fun, recycle, ...)
	if (test$nl < 1) stop("lapp does not like 'fun'")
	out <- rast(x, nlyr=test$nl)
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- lapply(1:length(x), function(s) readValues(x[s], b$row[i], b$nrows[i], 1, ncx, mat=TRUE))
		if (recycle) {
			v <- lapply(v, as.vector)
		}
		v <- do.call(fun, c(v, list(...)))
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	readStop(x)
	out <- writeStop(out)
	return(out)
}
)

