# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


setMethod("focal", signature(x="SpatRaster"), 
function(x, w=3, fun="sum", na.rm=TRUE, na.only=FALSE, fillvalue=NA, expand=FALSE, filename="",  ...)  {

	if (na.only && (!is.matrix(w))) {
		if (!na.rm) {
			warn("focal", "na.rm set to TRUE because na.only is TRUE")
			na.rm <- TRUE
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
	cpp <- FALSE
	txtfun <- .makeTextFun(fun)
	if (is.character(txtfun)) { 
		cpp <- TRUE
	}

	if (cpp) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$focal3(w, m, fillvalue, na.rm[1], na.only[1], txtfun, expand, opt)
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
		test <- try( apply(rbind(1:prod(w)), 1, fun, na.rm=FALSE), silent=TRUE )
		if (inherits(test, "try-error")) {
			test <- try( apply(rbind(1:prod(w)), 1, fun), silent=TRUE )
			if (!inherits(test, "try-error")) {
				usenarm = FALSE
			}
		}

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
		b <- writeStart(out, filename, n=msz*4, ...)
		
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
				if (usenarm) {
					v <- apply(v, 1, fun, na.rm=na.rm)
				} else {
					v <- apply(v, 1, fun)			
				}
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





setMethod("focalReg", signature(x="SpatRaster"), 
function(x, w=3, intercept=FALSE, na.rm=TRUE, fillvalue=NA, expand=FALSE, filename="",  ...)  {

	slope <- function(x, y) {
	  v <- cbind(x, y)
	  X <- cbind(1, v[,1])
	  invXtX <- solve(t(X) %*% X) %*% t(X)
	  (invXtX %*% v[,2])[2]
	}

	slopenarm <- function(x, y) {
	  v <- na.omit(cbind(x, y))
	  X <- cbind(1, v[,1])
	  invXtX <- solve(t(X) %*% X) %*% t(X)
	  (invXtX %*% v[,2])[2]
	}


	intslope <- function(x, y) {
	  v <- cbind(x, y)
	  X <- cbind(1, v[,1])
	  invXtX <- solve(t(X) %*% X) %*% t(X)
	  (invXtX %*% v[,2])
	}

	intslopenarm <- function(x, y) {
	  v <- na.omit(cbind(x, y))
	  X <- cbind(1, v[,1])
	  invXtX <- solve(t(X) %*% X) %*% t(X)
	  (invXtX %*% v[,2])
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

