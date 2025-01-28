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
		msg <- "cannot use 'fun'"
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {
			if (is.null(dim(vtst))) {
				msg <- paste0("cannot use 'fun'. The number of values returned is not divisible by the number of input cells (returning: ", length(vtst), ", expecting :", nr, ")")
			} else {
				msg <- paste0("cannot use 'fun'. The number of rows returned is not divisible by the number of input cells (returning: ", nrow(vtst), ", expecting: ", nr, ")")
			}
			nl <- -1
		}
	} else {
		if (is.null(dim(vtst))) {
			msg <- paste0("cannot use 'fun'. The number of values returned is less than the number of input cells.\n(returning: ", length(vtst), ", expecting: ", nr, ")\nPerhaps the function is not properly vectorized")
		} else {
			msg <- paste("cannot use 'fun'. The number of rows returned is less than the number of input cells.\n(returning:", nrow(vtst), ", expecting:", nr, ")\nPerhaps the function is not properly vectorized")
		}
		nl <- -1
	}
	if (nl < 0) {
		error("lapp", msg)
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

#	if (usenames) {
#		fnames <- names(formals(fun))
#		i <- names(x) %in% fnames
#		if (!all(i)) {
#			warn("lapp", paste(paste(names(x)[!i], collapse=","), "are not used", collapse=", "))
#			x <- x[[names(x) %in% fnames]]
#		}
#	}


	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		doclust <- TRUE
	} else if (cores > 1) {
		doclust <- TRUE
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores), add=TRUE)
	}

	readStart(x)
	on.exit(readStop(x), add=TRUE)
	ncx <- ncol(x)
	v <- readValues(x, round(0.51*nrow(x)), 1, 1, ncx, dataframe=TRUE)
	test <- .lapp_test(v, fun, usenames, ...)
	if (test$nl < 1) error("lapp", "cannot use 'fun'. The number of values returned is not divisible by the number of input cells")
	out <- rast(x, nlyrs=test$nl)
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	b <- writeStart(out, filename, overwrite, sources=sources(x), wopt=wopt)
	expected <- test$nl * ncx

	if (doclust) {
		ncores <- length(cores)
		export_args(cores, ..., caller="lapp")		
		cfun <- function(i, ...)  do.call(fun, i, ...)
		parallel::clusterExport(cores, "cfun", environment())
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, dataframe=TRUE)
			if (!usenames) colnames(v) <- NULL
			v <- split(v, rep(1:ncores, each=ceiling(nrow(v) / ncores))[1:nrow(v)])
			v <- unlist(parallel::parLapply(cores, v, cfun, ...))
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




.lapp_test_stack_call <- function(v, fun, recycle, ...) {
# figure out the shape of the output
	nms <- msg <- ""
	nr <- nrow(v[[1]])
	if (recycle) {
		v <- lapply(v, as.vector)
	}
	vtst <- try(do.call(fun, c(v, list(...))), silent=FALSE)
#	vtst2 <- try(apply(v, fun, ...), silent=TRUE)

	if (inherits(vtst, "try-error")) {
		nl <- -1
		msg <- "cannot use 'fun'"
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {
			if (is.null(dim(vtst))) {
				msg <- paste0("cannot use 'fun'. The number of values returned is not divisible by the number of input cells (returning: ", length(vtst), ", expecting :", nr, ")")
			} else {
				msg <- paste0("cannot use 'fun'. The number of rows returned is not divisible by the number of input cells (returning: ", nrow(vtst), ", expecting: ", nr, ")")
			}
			nl <- -1
		}
	} else {
		if (is.null(dim(vtst))) {
			msg <- paste0("cannot use 'fun'. The number of values returned is less than the number of input cells.\n(returning: ", length(vtst), ", expecting: ", nr, ")\nPerhaps the function is not properly vectorized.")
		} else {
			msg <- paste("cannot use 'fun'. The number of rows returned is less than the number of input cells.\n(returning:", nrow(vtst), ", expecting:", nr, ")\nPerhaps the function is not properly vectorized.")
		}
		nl <- -1
	}
	if (nl > 0) {
		if (is.matrix(vtst)) {
			nms <- colnames(vtst)
		}
	}
	list(nl=nl, names=nms, msg=msg)
}



.lapp_test_stack_mapp <- function(v, fun, recycle, ...) {
# figure out the shape of the output
	nms <- msg <- ""
	nr <- nrow(v[[1]])
	if (recycle) {
		v <- lapply(v, as.vector)
	}
	v <- lapply(v, function(i) data.frame(t(i)))
	vtst <- try(do.call(mapply, c(v, list(...), FUN=fun)), silent=FALSE)
	if (inherits(vtst, "try-error")) {
		return(list(nl=-10, names="", msg="cannot use 'fun'", trans=FALSE))
	}
	trans <- FALSE
	if (!is.null(dim(vtst))) {
		trans <- TRUE
		vtst <- as.vector(t(vtst))
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {
			if (is.null(dim(vtst))) {
				msg <- paste0("cannot use 'fun'. The number of values returned is not divisible by the number of input cells (returning: ", length(vtst), ", expecting :", nr, ")")
			} else {
				msg <- paste0("cannot use 'fun'. The number of rows returned is not divisible by the number of input cells (returning: ", nrow(vtst), ", expecting: ", nr, ")")
			}
			nl <- -1
		}
	} else {
		if (is.null(dim(vtst))) {
			msg <- paste0("cannot use 'fun'. The number of values returned is less than the number of input cells.\n(returning: ", length(vtst), ", expecting: ", nr, ")\nPerhaps the function is not properly vectorized.")
		} else {
			msg <- paste("cannot use 'fun'. The number of rows returned is less than the number of input cells.\n(returning:", nrow(vtst), ", expecting:", nr, ")\nPerhaps the function is not properly vectorized.")
		}
		nl <- -10
	}
	if (nl > 0) {
		if (is.matrix(vtst)) {
			nms <- colnames(vtst)
		}
	}
	list(nl=nl, names=nms, msg=msg, trans=trans)
}


setMethod("lapp", signature(x="SpatRasterDataset"),
function(x, fun, ..., usenames=FALSE, recycle=FALSE, cores=1, filename="", overwrite=FALSE, wopt=list())  {

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

	nms <- names(x)
	v <- lapply(1:length(x), function(i) readValues(x[i], round(0.51*nrx), 1, 1, ncx, mat=TRUE))
	if (usenames) names(v) <- nms
	mapp <- FALSE
	trans <- FALSE
	test <- .lapp_test_stack_call(v, fun, recycle, ...)
	if (test$nl < 1) {
		oldtst <- test
		test <- .lapp_test_stack_mapp(v, fun, recycle, ...)
		if (test$nl == 0) {
			error("lapp", paste0(oldtst$msg, "\n", test$msg))
		}
		mapp <- TRUE
		trans <- test$trans
	}
	
	out <- rast(x[1])
	nlyr(out) <- test$nl
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	nltot <- sum(nlyr(x)) + nlyr(out)
	fact <- max(4, 4 * nltot / nlyr(out))
	b <- writeStart(out, filename, overwrite, sources=unlist(sources(x)), wopt=wopt, n=fact)

	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		doclust <- TRUE
	} else if (cores > 1) {
		doclust <- TRUE
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores), add=TRUE)
	}
	if (mapp) {
		if (doclust) {
			warn("lapp", "no parallel method for this case")
		}
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, mat=TRUE)
			if (recycle) {
				v <- lapply(v, as.vector)
			}
			if (usenames) {
				names(v) <- nms
			}
			v <- lapply(v, function(j) data.frame(t(j)))
			v <- do.call(mapply, c(v, list(...), FUN=fun))
			if (test$trans) {
				v <- as.vector(t(v))
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		} 
	} else {
		if (doclust) {
			ncores <- length(cores)
			export_args(cores, ..., caller="lapp")		
			cfun <- function(i, ...)  do.call(fun, i, ...)
			parallel::clusterExport(cores, "cfun", environment())

			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, mat=TRUE)
				if (usenames) { names(v) <- nms } else  { names(v) <- NULL }
				
				if (recycle) {
					repl <- max(sapply(v, function(j) prod(dim(j))))
					v <- lapply(v, function(j) rep_len(as.vector(j), repl))
					j <- split(1:length(v[[1]]), rep(1:ncores, each=ceiling(length(v[[1]]) / ncores))[1:length(v[[1]])])
					v <- parallel::parLapply(cores, 1:ncores, function(x, ...) cfun(lapply(v, function(d) d[j[[x]]]), ...))
				} else {
					j <- split(1:nrow(v[[1]]), rep(1:ncores, each=ceiling(nrow(v[[1]]) / ncores))[1:nrow(v[[1]])])
					v <- parallel::parLapply(cores, 1:ncores, function(x, ...) cfun(lapply(v, function(d) d[j[[x]], , drop=FALSE]), ...))
				}

				v <- unlist(v)
				if (length(v) != (b$nrows[i] * test$nl * ncx)) {
					out <- writeStop(out)
					error("lapp", "output length of fun is not correct")
				}
				writeValues(out, v, b$row[i], b$nrows[i])
			}
		} else {
		
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, mat=TRUE)

				if (recycle) {
					v <- lapply(v, as.vector)
				}
				if (usenames) {
					names(v) <- nms
				}
				v <- do.call(fun, c(v, list(...)))
				writeValues(out, v, b$row[i], b$nrows[i])
			}
		}
	}
	writeStop(out)
}
)

