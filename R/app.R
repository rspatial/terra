

.cpp_funs <- c("sum", "mean", "median", "modal", "which", "which.min", "which.max", "min", "max", "prod", "any", "all", "sd", "std", "first")

setMethod("sapp", signature(x="SpatRaster"),
function(x, fun, ..., filename="", overwrite=FALSE, wopt=list())  {
	x <- lapply(x, function(r) fun(r, ...))
	x <- rast(x)
	if (filename[1] != "") {
		writeRaster(x, filename, overwrite, wopt=wopt)
	} else {
		tighten(x)
	}
}
)



setMethod("sapp", signature(x="SpatRasterDataset"),
function(x, fun, ..., filename="", overwrite=FALSE, wopt=list())  {
	x <- lapply(as.list(x), function(r) app(r, fun, ...))
	x <- rast(x)
	if (filename != "") {
		writeRaster(x, filename, overwrite, wopt=wopt)
	} else {
		tighten(x)
	}
}
)


export_args <- function(cores, ..., caller="app") {
	vals <- list(...)
	if (length(vals) < 1) return(NULL)
	nms <- names(vals)
	if (any(nms == "")) {
		error(caller, "additional arguments must be named when using multiple cores")
	}
	for (i in seq_along(vals)) {
		assign(nms[i], force(vals[i])) 
	}
	parallel::clusterExport(cores, nms, envir=environment())
}

## similar from predict
#	dots <- list(...)
#	if (length(dots) > 0) {
#		nms <- names(dots)
#		dotsenv <- new.env()
#		lapply(1:length(dots), function(i) assign(nms[i], dots[[i]], envir=dotsenv))
#		parallel::clusterExport(cls, nms, dotsenv)
#	}


setMethod("app", signature(x="SpatRaster"),
function(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())  {

	txtfun <- .makeTextFun(fun)
	if (inherits(txtfun, "character")) {
		if (txtfun %in% .cpp_funs) {
			opt <- spatOptions(filename, overwrite, wopt=wopt)
			na.rm <- isTRUE(list(...)$na.rm)
			x@ptr <- x@ptr$summary(txtfun, na.rm, opt)
			return(messages(x, "app"))
		}
	}
	fun <- match.fun(fun)
	out <- rast(x)
	nlyr(out) <- 1
	nc <- ncol(x)
	readStart(x)
	on.exit(readStop(x))
	nl <- nlyr(x)

	dots <- list(...)
	if (length(dots) > 0) {
		test <- any(sapply(dots, function(i) inherits(i, "SpatRaster")))
		if (test) {
			error("app", "additional arguments cannot be a SpatRaster")
		}
	}
# figure out the shape of the output by testing with up to 13 cells
	teststart <- max(1, 0.5 * nc - 6)
	testend <- min(teststart + 12, nc)
	ntest <- 1 + testend - teststart

	v <- readValues(x, round(0.51*nrow(x)), 1, teststart, ntest, mat=TRUE)
	usefun <- FALSE
	if (nl==1) {
		r <- fun(v, ...)
		usefun <- TRUE
	} else {
		r <- try(apply(v, 1, fun, ...), silent=TRUE)
		if (inherits(r, "try-error")) {
			rr <- try(fun(v, ...), silent=TRUE)
			if (inherits(rr, "try-error")) {
				error("app", paste0("cannot use this function\n", attr(r, "condition")))
			} else {
				usefun <- TRUE
				r <- rr
			}
		}
	}
	if (is.list(r)) {
		if (length(unique(sapply(r, length))) >  1) {
			error("app", "'fun' returns a list (should be numeric or matrix).\nPerhaps because returned values have different lengths due to NAs in input?")
		} else {
			error("app", "'fun' returns a list (should be numeric or matrix)")
		}
	}
	trans <- FALSE
	if (NCOL(r) > 1) {
		#? if ((ncol(r) %% ntest) == 0) {
		if (ncol(r) == ntest) {
			nlyr(out) <- nrow(r)
			trans <- TRUE
			nms <- rownames(r)
		} else if (nrow(r) == ntest) {
			nlyr(out) <- ncol(r)
			nms <- colnames(r)
		} else {
			error("app", "the number of values returned by 'fun' is not appropriate\n(it should be the product of the number of cells and and a positive integer)")
		}
		if (is.null(wopt$names)) {
			wopt$names <- nms
		}
	} else {
		if ((length(r) %% ntest) != 0) {
			error("app", "the number of values returned by 'fun' is not appropriate")
		} else {
			nlyr(out) <- length(r) / ntest
		}
	}

	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		doclust <- TRUE
	} else if (cores > 1) {
		doclust <- TRUE
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores), add=TRUE)
	}

	ncops <- nlyr(x) / nlyr(out)
	ncops <- ifelse(ncops > 1, ceiling(ncops), 1) * 4
	b <- writeStart(out, filename, overwrite, wopt=wopt, n=ncops, sources=sources(x))

	if (doclust) {
		ncores <- length(cores)
		export_args(cores, ...)
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
			icsz <- max(min(100, ceiling(b$nrows[i] / ncores)), b$nrows[i])
			r <- parallel::parRapply(cores, v, fun, ..., chunk.size=icsz)
			if (nlyr(out) > 1) {
				r <- matrix(r, ncol=nlyr(out), byrow=TRUE)
			}
			writeValues(out, r, b$row[i], b$nrows[i])
		}
	} else {
		if (usefun) {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
				r <- fun(v, ...)
				if (trans) {
					r <- t(r)
				}
				writeValues(out, r, b$row[i], b$nrows[i])
			}
		} else {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
				r <- apply(v, 1, fun, ...)
				if (trans) {
					r <- t(r)
					#r <- as.vector(r)
				}
				writeValues(out, r, b$row[i], b$nrows[i])
			}
		}
	}
	writeStop(out)
}
)



.app_test_stack <- function(v, fun, ncols, ...) {
# figure out the shape of the output
	nms = ""
	nr <- nrow(v[[1]])
	v <- lapply(v, as.vector)
	v <- do.call(cbind, v)
	r <- apply(v, 1, fun, ...)
	if (inherits(r, "try-error")) {
		nl <- -1
	}

	trans <- FALSE
	if (NCOL(r) > 1) {
		#? if ((ncol(r) %% nc) == 0) {
		if (ncol(r) == ncols) {
			nl <- nrow(r)
			trans <- TRUE
		} else if (nrow(r) == ncols) {
			nl <- ncol(r)
		} else {
			error("app", "cannot handle 'fun'")
		}
	} else if (length(r) >= nr) {
		if ((length(r) %% nr) == 0) {
			nl <- length(r) / nr
		} else {
			nl <- -1
		}
	} else {
		nl <- -1
	}
	if (is.matrix(r)) {
		nms <- colnames(r)
	}
	list(nl=nl, trans=trans, names=nms)
}



.app_test_stack <- function(v, fun, ncols, ...) {
# figure out the shape of the output
	nms = ""
	nr <- nrow(v[[1]])
	v <- lapply(v, as.vector)
	v <- do.call(cbind, v)
	r <- apply(v, 1, fun, ...)
	if (inherits(r, "try-error")) {
		nl <- -1
	}

	trans <- FALSE
	if (NCOL(r) > 1) {
		#? if ((ncol(r) %% nc) == 0) {
		if (ncol(r) == ncols) {
			nl <- nrow(r)
			trans <- TRUE
		} else if (nrow(r) == ncols) {
			nl <- ncol(r)
		} else {
			error("app", "'fun' is not appropriate")
		}
	} else if (length(r) >= nr) {
		if ((length(r) %% nr) == 0) {
			nl <- length(r) / nr
		} else {
			nl <- -1
		}
	} else {
		nl <- -1
	}
	if (is.matrix(r)) {
		nms <- colnames(r)
	}
	list(nl=nl, trans=trans, names=nms)
}



setMethod("app", signature(x="SpatRasterDataset"),
function(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())  {

	txtfun <- .makeTextFun(match.fun(fun))
	if (inherits(txtfun, "character")) {
		if (txtfun %in% .cpp_funs) {
			opt <- spatOptions(filename, overwrite, wopt=wopt)
			narm <- isTRUE(list(...)$na.rm)
			r <- rast()
			r@ptr <- x@ptr$summary(txtfun, narm, opt)
			return (messages(r, "app") )
		}
	}

	if (missing(fun)) error("app", "'fun' is missing")

	ncx <- ncol(x[1])
	nrx <- nrow(x[1])
	readStart(x)
	on.exit(readStop(x))

	v <- lapply(1:length(x), function(i) readValues(x[i], round(0.51*nrx), 1, 1, ncx, mat=TRUE))
	test <- .app_test_stack(v, fun, ncx, ...)
	if (test$nl < 1) error("app", "cannot find 'fun'")
	out <- rast(x[1], nlyrs=test$nl)
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}

	nc <- (nlyr(x[1]) * length(x)) / nlyr(out)
	nc <- ifelse(nc > 1, ceiling(nc), 1) * 3
	b <- writeStart(out, filename, overwrite, wopt=wopt, n=nc, sources=unlist(sources(x)))


	if (inherits(cores, "cluster")) {
		doclust <- TRUE
	} else if (cores > 1) {
		doclust <- TRUE
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores), add=TRUE)
	} else {
		doclust <- FALSE
	}

	if (doclust) {
		ncores <- length(cores)
		export_args(cores, ...)
		for (i in 1:b$n) {
			v <- lapply(1:length(x), function(s) as.vector(readValues(x[s], b$row[i], b$nrows[i], 1, ncx, mat=TRUE)))
			v <- do.call(cbind, v)
			icsz <- max(min(100, ceiling(b$nrows[i] / ncores)), b$nrows[i])
			r <- parallel::parRapply(cores, v, fun, ..., chunk.size=icsz)
			if (test$trans) {
				r <- t(r)
			}
			writeValues(out, r, b$row[i], b$nrows[i])

		}
	} else {
		for (i in 1:b$n) {
			v <- lapply(1:length(x), function(s) as.vector(readValues(x[s], b$row[i], b$nrows[i], 1, ncx, mat=TRUE)))
			r <- apply(do.call(cbind, v), 1, fun, ...)
			if (test$trans) {
				r <- t(r)
			}
			writeValues(out, r, b$row[i], b$nrows[i])
		}
	}
	#readStop(x)
	writeStop(out)
}
)


