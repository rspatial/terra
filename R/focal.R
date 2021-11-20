# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, fun="sum", ..., na.only=FALSE, fillvalue=NA, expand=FALSE, filename="", overwrite=FALSE, wopt=list())  {


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
	cpp <- FALSE
	txtfun <- .makeTextFun(fun)
	if (is.character(txtfun)) { 
		cpp <- TRUE
	}

	if (cpp) {
		opt <- spatOptions(filename, overwrite, wopt)
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
		
		usenarm = TRUE
		test <- apply(rbind(1:prod(w)), 1, fun, ...)

		nl <- nlyr(x)
		outnl <- nl * length(test)
		transp <- FALSE
		nms <- NULL
		if (isTRUE(nrow(test) > 1)) {
			transp <- TRUE
			nms <- rownames(test)
		} else if (isTRUE(ncol(test) > 1)) {
			nms <- rownames(test)
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
function(x, w=3, intercept=FALSE, na.rm=TRUE, fillvalue=NA, expand=FALSE, filename="",  ...)  {

	slope <- function(x, y) {
		v <- cbind(x, y)
		X <- cbind(1, v[,1])
		XtX <- t(X) %*% X
		if (det(XtX) == 0) {
			return(NA)
		}
		invXtX <- solve(XtX) %*% t(X)
		(invXtX %*% v[,2])[2]
	}

	slopenarm <- function(x, y) {
		v <- na.omit(cbind(x, y))
		if (nrow(v) == 0) {
			return(NA)
		}
		X <- cbind(1, v[,1])
		XtX <- t(X) %*% X
		if (det(XtX) == 0) {
			return(NA)
		}
		invXtX <- solve(XtX) %*% t(X)
		(invXtX %*% v[,2])[2]
	}


	intslope <- function(x, y) {
		v <- cbind(x, y)
		X <- cbind(1, v[,1])
		XtX <- t(X) %*% X
		if (det(XtX) == 0) {
			return(NA)
		}
		invXtX <- solve(XtX) %*% t(X)
		invXtX %*% v[,2]
	}

	intslopenarm <- function(x, y) {
		v <- na.omit(cbind(x, y))
		if (nrow(v) == 0) {
			return(NA)
		}
		X <- cbind(1, v[,1])
		XtX <- t(X) %*% X
		if (det(XtX) == 0) {
			return(NA)
		}
		invXtX <- solve(XtX) %*% t(X)
		invXtX %*% v[,2]
	}

	if (nlyr(x) != 2) error("focalReg", "x must have 2 layers")
	outnl <- 1 + isTRUE(intercept)	
	if (na.rm) {
		if (outnl == 2) {
			fun = function(x, y) try(intslopenarm(x, y), silent=TRUE)
		} else {
			fun = function(x, y) try(slopenarm(x, y), silent=TRUE)		
		}
	} else {
		if (outnl == 2) {
			fun = function(x, y) try(intslope(x, y), silent=TRUE)
		} else {
			fun = function(x, y) try(slope(x, y), silent=TRUE)		
		}
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
	dow <- !isTRUE(all(m == 1))
	if (any(is.na(m))) {
		k <- !is.na(m)
		mm <- m[k]
		msz <- sum(k)
	}
	out <- rast(x, nlyr=outnl)
	if (outnl == 2) {
		names(out) <- c("intercept", "slope")
	} else {
		names(out) <- "slope"
	}
	b <- writeStart(out, filename, n=msz*4, ...)

	for (i in 1:b$n) {
		X <- focalValues(x[[1]], w)
		Y <- focalValues(x[[2]], w)
		if (dow) {
			if (any(is.na(m))) {
				X <- X[k] * mm
				Y <- Y[k] * mm
			} else {
				X <- X * m
				Y <- Y * m
			}
		}
		v <- sapply(1:nrow(X), function(i) fun(X[i,], Y[i,]))
		if (outnl == 2) {
			v <- as.vector(t(v))
		}
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)

