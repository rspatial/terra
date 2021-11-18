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
function(x, fun, ..., usenames=FALSE, cores=1, filename="", overwrite=FALSE, wopt=list())  {

	fun <- match.fun(fun)
	dots <- list(...)
	if (any(sapply(dots, function(i) inherits(i, "SpatRaster")))) {
		error("lapp", "only 'x' can be a SpatRaster")
		# otherwise .lapp_test may crash! 
	}

	if (usenames) {
		fnames <- names(formals(fun))
		x <- x[[names(x) %in% fnames]]
	}


	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		doclust <- TRUE
		ncores <- length(cores)		
	} else if (cores > 1) {
		doclust <- TRUE
		ncores <- cores
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores), add=TRUE)	
	}

	readStart(x)
	on.exit(readStop(x), add=TRUE)
	ncx <- ncol(x)
	v <- readValues(x, round(0.51*nrow(x)), 1, 1, ncx, dataframe=TRUE)
	test <- .lapp_test(v, fun, usenames, ...)
	if (test$nl < 1) error("lapp", "I do not like 'fun' :(")
	out <- rast(x, nlyrs=test$nl)
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	b <- writeStart(out, filename, overwrite, wopt=wopt)
	expected <- test$nl * ncx

	if (doclust) {
		cfun <- \(i, ...)  do.call(fun, i, ...)
		parallel::clusterExport(cores, "cfun", environment())
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, dataframe=TRUE)
			if (!usenames) colnames(v) <- NULL
			v <- split(v, rep(1:ncores, each=ceiling(nrow(v) / ncores))[1:nrow(v)])
			v <- unlist(parallel::parLapply(cores, v, cfun))
			if (length(v) != (expected * b$nrows[i])) {
				out <- writeStop(out)
				error("lapp", "output length of fun is not correct")
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	} else {
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, dataframe=TRUE)
			if (!usenames) colnames(v) <- NULL
			v <- do.call(fun, c(v, list(...)))
			if (length(v) != (expected * b$nrows[i])) {
				out <- writeStop(out)
				error("lapp", "output length of fun is not correct")
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		}
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



setMethod("lapp", signature(x="SpatRasterDataset"), 
function(x, fun, ..., recycle=FALSE, filename="", overwrite=FALSE, wopt=list())  {

	fun <- match.fun(fun)
	dots <- list(...)
	if (any(sapply(dots, function(i) inherits(i, "SpatRasterDataset")))) {
		error("lapp", "only 'x' can be a SpatRasterDataset")
		# otherwise .lapp_test_stack fails
	}

	ncx <- ncol(x[1])
	nrx <- nrow(x[1])
	readStart(x)
	on.exit(readStop(x))

	v <- lapply(1:length(x), function(i) readValues(x[i], round(0.51*nrx), 1, 1, ncx, mat=TRUE))
	test <- .lapp_test_stack(v, fun, recycle, ...)
	if (test$nl < 1) error("lapp", "cannot use 'fun'")
	out <- rast(x)
	nlyr(out) <- test$nl
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	nltot <- sum(nlyr(x)) + nlyr(out)
	fact <- max(4, 4 * nltot / nlyr(out))
	b <- writeStart(out, filename, overwrite, wopt=wopt, n=fact)
	for (i in 1:b$n) {
		v <- lapply(1:length(x), function(s) readValues(x[s], b$row[i], b$nrows[i], 1, ncx, mat=TRUE))
		if (recycle) {
			v <- lapply(v, as.vector)
		}
		v <- do.call(fun, c(v, list(...)))
		writeValues(out, v, b$row[i], b$nrows[i])
	}
#	readStop(x)
	out <- writeStop(out)
	return(out)
}
)

