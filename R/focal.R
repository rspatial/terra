# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, fun="sum", ..., na.only=FALSE, fillvalue=NA, expand=FALSE, filename="", overwrite=FALSE, wopt=list())  {

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
		x@ptr <- x@ptr$focal3(w, m, fillvalue, narm, na.only[1], txtfun, expand, opt)
		messages(x, "focal")
		return(x)

	} else {
		msz <- prod(w)
		readStart(x)
		on.exit(readStop(x))
		dow <- !isTRUE(all(m == 1))
		if (any(is.na(m))) {
			k <- !is.na(m)
			mm <- m[k]
			msz <- sum(k)
		}
		
		usenarm <- TRUE
		test <- apply(rbind(1:prod(w)), 1, fun, ...)

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
		
		for (i in 1:b$n) {
			vv <- NULL
			for (j in 1:nl) {
				if (nl > 1) {
					v <- x[[j]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i])
				} else {
					v <- x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i])
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

				if (na.only) {
					vv <- readValues(x, b$row[i], b$nrows[i])
					j <- !is.na(vv)
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


setMethod("focalCpp", signature(x="SpatRaster"), 
function(x, w=3, fun, ..., fillvalue=NA, expand=FALSE, filename="", overwrite=FALSE, wopt=list())  {

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

	test <- fun(1:msz, ..., ni=1, nw=msz)
	nl <- nlyr(x)
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
				v <- x[[j]]@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i])
			} else {
				v <- x@ptr$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i])
			}
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
	out <- writeStop(out)
	return(out)
}
)




setMethod("focalReg", signature(x="SpatRaster"), 
function(x, w=3, na.rm=TRUE, fillvalue=NA, expand=FALSE, filename="",  ...)  {

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
function(x, w=3, fun, ..., fillvalue=NA, expand=FALSE, filename="", overwrite=FALSE, wopt=list()) {

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

