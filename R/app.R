
setMethod("app", signature(x="SpatRaster"), 
function(x, fun, ..., nodes=1, filename="", overwrite=FALSE, wopt=list())  {

	txtfun <- .makeTextFun(match.fun(fun))
	if (inherits(txtfun, "character")) { 
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
	if (!readStart(x)) { stop(x@ptr$messages$getError()) }
	on.exit(readStop(x))

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
			nms <- rownames(r)
		} else if (nrow(r) == nc) {
			nlyr(out) <- ncol(r)
			nms <- colnames(r)
		} else {
			stop("cannot handle this function")
		}
		if (is.null(wopt$names)) {
			wopt$names <- nms
		}
	}

	b <- writeStart(out, filename, overwrite, wopt)
	if (nodes > 1) {
		cls <- parallel::makeCluster(nodes)
		on.exit(parallel::stopCluster(cls))
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
			icsz <- max(min(100, ceiling(b$nrows[i] / nodes)), b$nrows[i])
			r <- parallel::parRapply(cls, v, fun, ..., chunk.size=icsz)
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
			stop("cannot handle this function")
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



setMethod("app", signature(x="SpatDataSet"), 
function(x, fun, ..., nodes=1, filename="", overwrite=FALSE, wopt=list())  {

	txtfun <- .makeTextFun(match.fun(fun))
	if (inherits(txtfun, "character")) { 
		if (txtfun %in% c("max", "min", "mean", "range", "prod", "sum", "any", "all")) {
			opt <- .runOptions(filename, overwrite, wopt)
			narm <- isTRUE(list(...)$na.rm)	
			r <- rast()
			r@ptr <- x@ptr$summary(txtfun, narm, .terra_environment$options@ptr)
			return (show_messages(r, "app") )
		}		
	}

	stopifnot(!missing(fun))
	
	ncx <- ncol(x[1])
	nrx <- nrow(x[1])
	if (!readStart(x)) { stop(x@ptr$messages$getError()) }
	on.exit(readStop(x))

	v <- lapply(1:length(x), function(i) readValues(x[i], round(0.5*nrx), 1, 1, ncx, mat=TRUE))
	test <- .app_test_stack(v, fun, ncx, ...)
	if (test$nl < 1) stop("app is not having 'fun'")
	out <- rast(x[1], nlyr=test$nl)
	if (length(test$names == test$nl)) {
		if (is.null(wopt$names)) wopt$names <- test$names
	}
	b <- writeStart(out, filename, overwrite, wopt)

	if (nodes > 1) {
		cls <- parallel::makeCluster(nodes)
		on.exit(parallel::stopCluster(cls))
		for (i in 1:b$n) {
			v <- lapply(1:length(x), function(s) as.vector(readValues(x[s], b$row[i], b$nrows[i], 1, ncx, mat=TRUE)))
			v <- do.call(cbind, v)
			icsz <- max(min(100, ceiling(b$nrows[i] / nodes)), b$nrows[i])
			r <- parallel::parRapply(cls, v, fun, ..., chunk.size=icsz)
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
	readStop(x)
	out <- writeStop(out)
	return(out)
}
)



