# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"),
function(x, w=3, fun="sum", ..., na.policy="all", fillvalue=NA, expand=FALSE, silent=TRUE, filename="", overwrite=FALSE, wopt=list())  {

	na.only <- list(...)$na.only
	if (!is.null(na.only)) {
		warn("focal", "use 'na.policy' instead of 'na.only'")
		na.policy <- "only"
	}
	na.policy <- match.arg(tolower(na.policy), c("all", "only", "omit"))
	na.only <- na.policy == "only"
	na.omit <- na.policy == "omit"

	if (!is.numeric(w)) {
		error("focal", "w should be numeric vector or matrix")
	}
	txtfun <- .makeTextFun(fun)

	if (is.matrix(w)) {
		m <- as.vector(t(w))
		if (!all(m %in% c(0, 1, NA))) {
			if (isTRUE(list(...)$na.rm)) {
				if (txtfun != "sum") {
					error("focal", 'with "na.rm=TRUE" and weights other than 0, 1, or NA, only fun="sum" is allowed')
				}
			}
		}
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}
	cpp <- FALSE
	txtfun <- .makeTextFun(fun)
	if (is.character(txtfun)) {
		if (is.null(wopt$names)) {
			wopt$names <- paste0("focal_", txtfun)
		}
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		narm <- isTRUE(list(...)$na.rm)
		if (na.only && (!narm)) {
			error("focal", "combining 'na.only=TRUE' with 'na.rm=FALSE' has no effect")
		}
		x@ptr <- x@ptr$focal(w, m, fillvalue, narm, na.only, na.omit, txtfun, expand, opt)
		messages(x, "focal")
		return(x)

	} else {
		if (expand) {
			warn(focal, "expand is ignored for non-standard functions")
		}
		checkNA <- na.only || na.omit

		msz <- prod(w)
		dow <- !isTRUE(all(m == 1))
		if (any(is.na(m))) {
			k <- !is.na(m)
			mm <- m[k]
			msz <- sum(k)
		}

		v <- focalValues(x, w, trunc(nrow(x)/2), 1)[ncol(x)/2, ,drop=FALSE]
		if (dow) {
			if (any(is.na(m))) {
				v <- v[k] * mm
			} else {
				v <- v * m
			}
		}
		test <- try(apply(v, 1, fun, ...), silent=silent)
		if (inherits(test, "try-error")) {
			error("focal", "test failed")
		}

		readStart(x)
		on.exit(readStop(x))
		nl <- nlyr(x)
		outnl <- nl * length(test)
		transp <- FALSE
		nms <- NULL
		if (isTRUE(nrow(test) > 1)) {
			transp <- TRUE
			nms <- rownames(test)
		} else if (isTRUE(ncol(test) > 1)) {
			nms <- colnames(test)
		}

		out <- rast(x, nlyr=outnl)
		if (!is.null(nms)) {
			names(out) <- nms
		}
		b <- writeStart(out, filename, overwrite, n=msz*4, wopt=wopt)
		opt <- spatOptions()

		for (i in 1:b$n) {
			vv <- NULL
			for (j in 1:nl) {
				if (nl > 1) {
					v <- x[[j]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
				} else {
					v <- x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
				}
				if (dow) {
					if (any(is.na(m))) {
						v <- v[k] * mm
					} else {
						v <- v * m
					}
				}
				v <- matrix(v, ncol=msz, byrow=TRUE)
				v <- apply(v, 1, fun, ...)
				if (transp) {
					v <- t(v)
				}

				if (checkNA) {
					vv <- readValues(x, b$row[i], b$nrows[i])
					if (na.only) {
						j <- !is.na(vv)
					} else {
						j <- is.na(vv)
					}
					v[j] <- vv[j]
				}
				if (nl > 1) {
					if (outnl > 1) {
						vv <- rbind(vv, v)
					} else {
						vv <- c(vv, v)
					}
				}
			}
			#if (bip) {
			#	v <- matrix(as.vector(v), ncol=ncol(v), byrow=TRUE)
			#}
			if (nl > 1) {
				writeValues(out, vv, b$row[i], b$nrows[i])
			} else {
				writeValues(out, v, b$row[i], b$nrows[i])
			}
		}
		out <- writeStop(out)
		return(out)
	}
}
)

setMethod("focal3D", signature(x="SpatRaster"),
function(x, w=3, fun=mean, ..., na.policy="all", fillvalue=NA, pad=FALSE, padvalue=fillvalue, expand=FALSE, silent=TRUE, filename="", overwrite=FALSE, wopt=list()) {


	na.policy <- match.arg(tolower(na.policy), c("all", "only", "omit"))
	na.only <- na.policy == "only"
	na.omit <- na.policy == "omit"
	checkNA <- na.only || na.omit

	if (!(inherits(w, "numeric") || inherits(w, "array"))) {
		error("focal3D", "w should be numeric vector or array")
	}

	if (is.array(w)) {
		if (length(dim(w)) != 3) {
			error("focal3D", "the weights array must have three dimensions")
		}
		m <- as.vector(w)
		w <- dim(w)
	} else {
		w <- rep_len(w, 3)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}
	if (w[3] > nlyr(x)) {
		error("focal", "the third weights dimension is larger than nlyr(x)")
	}
	
	msz <- prod(w)
	dow <- !isTRUE(all(m == 1))
	rna <- FALSE
	if (any(is.na(m))) {
		rna <- TRUE
		kna <- !is.na(m)
		m <- m[kna]
		msz <- sum(kna)
	} 

	opt <- spatOptions()
	halfway <- floor(w[3]/2)
	v <- lapply(1:w[3], function(i) focalValues(x[[i]], w[1:2], trunc(nrow(x)/2), 1)[ncol(x)/2, ,drop=FALSE])
	v <- t(do.call(cbind, v))
	if (dow) {
		if (rna) {
			v <- v[kna,,drop=FALSE] * m
		} else {
			v <- v * m
		}
	}
	vout <- try(apply(v, 2, fun, ...), silent=silent)
	if (inherits(vout, "try-error")) {
		error("focal", "test failed")
	}

	readStart(x)
	on.exit(readStop(x))
	nl <- nlyr(x)
	transp <- FALSE
	nms <- NULL
	if (isTRUE(nrow(vout) > 1)) {
		transp <- TRUE
		nms <- rownames(vout)
	} else if (isTRUE(ncol(vout) > 1)) {
		nms <- colnames(vout)
	}

	if (pad || expand) {
		startlyr = 1;
		endlyr = nl;
		outnl <- nl * length(vout)
	} else {
		startlyr = halfway+1;
		endlyr = nl-halfway;	
		outnl <- (1+endlyr-startlyr) * length(vout)
	}

	out <- rast(x, nlyr=outnl)
	if (!is.null(nms)) {
		names(out) <- nms
	}
	b <- writeStart(out, filename, overwrite, n=msz*4, wopt=wopt)

	nread <- prod(w[1:2])

	for (i in 1:b$n) {
		nc <- b$nrows[i]*ncol(x)
		vv <- NULL
		if (expand) {
			v <- list(matrix(x[[1]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			v <- do.call(rbind, rep(v, halfway+1))
			for (k in 2:(1+halfway)) {
				v <- rbind(v,  matrix(x[[k]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			}
		} else if (pad) {
			v <- matrix(padvalue, ncol=b$nrows[i]*ncol(x), nrow=nread*halfway)
			for (k in 1:(1+halfway)) {
				v <- rbind(v,  matrix(x[[k]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			}
		} else {
			v <- lapply(1:w[3], 
				function(k) matrix(x[[k]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			v <- do.call(rbind, v)
		}
		for (j in startlyr:endlyr) {
			if (j > startlyr) {
				v <- v[-c(1:nread), ]
				k <- j + halfway
				if (k > nl) {
					if (pad) {
						v <- rbind(v, matrix(padvalue, nrow=nread, ncol=ncol(v)))
					} else {
						v <- rbind(v, v[(nrow(v)-nread):nrow(v), ])
					}
				} else {
					v <- rbind(v, matrix(x[[k]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=ncol(v)))
				}
			}
			if (dow) {
				if (rna) {
					vout <- apply(v[kna,] * m, 2, fun, ...)
				} else {
					vout <- apply(v * m, 2, fun, ...)
				}
			} else {
				vout <- apply(v, 2, fun,...)
			}
			
			if (transp) {
				vout <- t(vout)
			}
			if (checkNA) {
				vv <- readValues(x, b$row[i], b$nrows[i])
				if (na.only) {
					j <- !is.na(vv)
				} else {
					j <- is.na(vv)
				}
				vout[j] <- vv[j]
			}
			vv <- c(vv, as.vector(vout))
		}
		writeValues(out, vv, b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)



setMethod("focalCpp", signature(x="SpatRaster"),
function(x, w=3, fun, ..., fillvalue=NA, silent=TRUE, filename="", overwrite=FALSE, wopt=list())  {

	if (!(all(c("ni", "nw") %in% names(formals(fun))))) {
		error("focalRaw", 'fun must have an argument "ni"')
	}

	if (!is.numeric(w)) {
		error("focal", "w should be numeric vector or matrix")
	}
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}

	msz <- prod(w)
	readStart(x)
	on.exit(readStop(x))
	dow <- !isTRUE(all(m == 1))
	if (any(is.na(m))) {
		k <- !is.na(m)
		mm <- m[k]
		msz <- sum(k)
	}

	opt <- spatOptions()
	nl <- nlyr(x)
	v <- x@ptr$focalValues(w, fillvalue, max(0, trunc(nrow(x)/2)), 1, opt)[1:prod(w)]
	if (dow) {
		if (any(is.na(m))) {
			v <- v[k] * mm
		} else {
			v <- v * m
		}
	}
	test <- try(fun(v, ..., ni=1, nw=msz), silent=silent)
	if (inherits(test, "try-error")) {
		error("focalCpp", "test failed")
	}
	outnl <- nl * length(test)
	if (is.null(wopt$names )) {
		wopt$names <- colnames(test)
	}
	out <- rast(x, nlyr=outnl)
	b <- writeStart(out, filename, overwrite, n=msz*4, wopt=wopt)

	nc <- ncol(out)
	for (i in 1:b$n) {
		vv <- NULL
		for (j in 1:nl) {
			if (nl > 1) {
				v <- x[[j]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
			} else {
				v <- x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
			}
			sst <- messages(x)
			if (dow) {
				if (any(is.na(m))) {
					v <- v[k] * mm
				} else {
					v <- v * m
				}
			}
			v <- fun(v, ..., ni=b$nrows[i]*nc, nw=msz)
			if (nl > 1) {
				if (outnl > 1) {
					vv <- rbind(vv, v)
				} else {
					vv <- c(vv, v)
				}
			}
		}
		if (nl > 1) {
			writeValues(out, vv, b$row[i], b$nrows[i])
		} else {
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	}
	writeStop(out)
}
)




setMethod("focalReg", signature(x="SpatRaster"),
function(x, w=3, na.rm=TRUE, fillvalue=NA, filename="",  ...)  {

	ols <- function(x, y) {
		v <- cbind(y, x)
		X <- cbind(1, v[,-1])
		XtX <- t(X) %*% X
		if (det(XtX) == 0) {
			return(NA)
		}
		invXtX <- solve(XtX) %*% t(X)
		invXtX %*% v[,1]
	}

	ols_narm <- function(x, y) {
		v <- na.omit(cbind(y, x))
		if (nrow(v) < (NCOL(x) + 1)) {
			return(NA)
		}
		X <- cbind(1, v[,-1])
		XtX <- t(X) %*% X
		if (det(XtX) == 0) {
			return(NA)
		}
		invXtX <- solve(XtX) %*% t(X)
		invXtX %*% v[,1]
	}

	nl <- nlyr(x)
	if (nl < 2) error("focalReg", "x must have at least 2 layers")

	if (!is.numeric(w)) {
		error("focalReg", "w should be numeric vector or matrix")
	}
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}
	msz <- prod(w)
	dow <- !isTRUE(all(m == 1))
	isnam <- FALSE
	if (any(is.na(m))) {
		k <- !is.na(m)
		mm <- m[k]
		msz <- sum(k)
		isnam <- TRUE
	}
	out <- rast(x)

	if (na.rm) {
		fun = function(x, y) {
			x <- try(ols_narm(x, y), silent=TRUE)
		}
		#fun = ols_narm
	} else {
		fun = function(x, y) try(ols(x, y), silent=TRUE)
		#fun = ols
	}
	names(out) <- paste0("B", 0:(nl-1))
	b <- writeStart(out, filename, n=msz*4, ...)
	ry <- x[[1]]
	rx <- x[[-1]]

	if (nl == 2) {
		for (i in 1:b$n) {
			Y <- focalValues(ry, w, b$row[i], b$nrows[i], fillvalue)
			X <- focalValues(rx, w, b$row[i], b$nrows[i], fillvalue)
			if (dow) {
				if (isnam) {
					Y <- Y[k] * mm
					X <- X[k] * mm
				} else {
					Y <- Y * m
					X <- X * m
				}
			}
			v <- t(sapply(1:nrow(Y), function(i) fun(X[i,], Y[i,])))
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	} else {

		for (i in 1:b$n) {
			Y <- focalValues(ry, w, b$row[i], b$nrows[i], fillvalue)
			if (dow) {
				if (isnam) {
					Y <- Y[k] * mm
				} else {
					Y <- Y * m
				}
			}
			X <- list()
			for (j in 1:(nl-1)) {
				X[[j]] <- focalValues(rx[[j]], w, b$row[i], b$nrows[i], fillvalue)
				if (dow) {
					if (any(is.na(m))) {
						X[[j]] <- X[[j]][k] * mm
					} else {
						X[[j]] <- X[[j]] * m
					}
				}
			}
			v <- list()
			for (p in 1:nrow(Y)) {
				xlst <- list()
				for (j in 1:(nl-1)) {
					xlst[[j]] <- X[[j]][p,]
				}
				pX <- do.call(cbind, xlst)
				v[[p]] <- fun(pX, Y[p,])
			}
			v <- t(do.call(cbind, v))
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	}

	out <- writeStop(out)
	return(out)
}
)



setMethod("focalCor", signature(x="SpatRaster"),
function(x, w=3, fun, ..., fillvalue=NA, filename="", overwrite=FALSE, wopt=list()) {

	nl <- nlyr(x)
	if (nl < 2) error("focalCor", "x must have at least 2 layers")

	if (!is.numeric(w)) {
		error("focalCor", "w should be numeric vector or matrix")
	}
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}
	msz <- prod(w)
	dow <- !isTRUE(all(m == 1))
	isnam <- FALSE
	if (any(is.na(m))) {
		k <- !is.na(m)
		mm <- m[k]
		msz <- sum(k)
		isnam <- TRUE
	}

	test <- do.call(fun, list(1:prod(w), prod(w):1), ...)
	if (is.null(wopt$names )) {
		wopt$names <- colnames(test)
	}

	outnl <- (nlyr(x) - 1) * length(test)
	out <- rast(x, nlyr=outnl)

	b <- writeStart(out, filename, n=msz*4, ...)

	v <- list()
	for (i in 1:b$n) {
		Y <- focalValues(x[[1]], w, b$row[i], b$nrows[i], fillvalue)
		if (dow) {
			if (isnam) {
				Y <- Y[k] * mm
			} else {
				Y <- Y * m
			}
		}
		for (j in 2:nlyr(x)) {
			X <- Y
			Y <- focalValues(x[[j]], w, b$row[i], b$nrows[i], fillvalue)
			if (dow) {
				if (isnam) {
					Y <- Y[k] * mm
				} else {
					Y <- Y * m
				}
			}
			v[[j-1]] <- t(sapply(1:nrow(Y), function(i, ...) fun(X[i,], Y[i,], ...)))
		}
		v <- do.call(cbind, v)
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)
