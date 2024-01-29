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
		#if (!all(m %in% c(0, 1, NA))) {
			#if (isTRUE(list(...)$na.rm)) {
			#	if (txtfun != "sum") {
					# error("focal", 'with "na.rm=TRUE" and weights other than 0, 1, or NA, only fun="sum" is allowed')
			#	}
			#}
		#}
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}
	cpp <- FALSE

	if (is.character(txtfun)) {
		if (is.null(wopt$names)) {
			wopt$names <- paste0("focal_", txtfun)
		}
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		if (na.only) {
			narm <- TRUE
		} else {
			narm <- isTRUE(list(...)$na.rm)
		}
		x@cpp <- x@cpp$focal(w, m, fillvalue, narm, na.only, na.omit, txtfun, expand, opt)
		messages(x, "focal")
		return(x)

	} else {
		if (expand) {
			warn("focal", "expand is ignored for functions that are not 'built-in'")
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
				v <- v[,k,drop=FALSE] * mm
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
		out <- rast(x, nlyr=outnl)

		transp <- FALSE
		nms <- NULL
		if (isTRUE(nrow(test) > 1)) {
			transp <- TRUE
			nms <- rownames(test)
			if (nl > 1) {
				nms <- paste0(rep(names(x), each=nl), "_", rep(nms, nl))
			}
		} else if (isTRUE(ncol(test) > 1)) {
			nms <- colnames(test)
			if (nl > 1) {
				nms <- paste0(rep(names(x), each=nl), "_", rep(nms, nl))
			}
		} else {
			nms <- names(x)
			if (nl > 1) {
				nms <- paste0(rep(names(x), each=nl), "_", rep(1:length(test), nl))
			}
		}
		if (length(nms) == nlyr(out)) {
			names(out) <- nms
		}
		b <- writeStart(out, filename, overwrite, n=msz+2, sources=sources(x), wopt=wopt)
		opt <- spatOptions()

		for (i in 1:b$n) {
			vv <- NULL
			for (j in 1:nl) {
				if (nl > 1) {
					v <- x[[j]]@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
				} else {
					v <- x@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
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
					if (nl > 1) {
						mv <- readValues(x[[j]], b$row[i], b$nrows[i])
					} else {
						mv <- readValues(x, b$row[i], b$nrows[i])
					}
					
					if (na.only) {
						k <- !is.na(mv)
					} else {
						k <- is.na(mv)
					}
					v[k] <- mv[k]
				}
				if (nl > 1) {
					if (outnl > 1) {
						vv <- cbind(vv, v)
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
		error("focal3D", "the third weights dimension is larger than nlyr(x)")
	}
	if (any((w %% 2) == 0)) {
		error("focal3D", "w must be odd sized in all dimensions")
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
		if (length(nms) <= outnl) {
			names(out) <- rep_len(nms, outnl) 
		}
	}
	b <- writeStart(out, filename, overwrite, n=msz+2, sources=sources(x), wopt=wopt)

	nread <- prod(w[1:2])

	for (i in 1:b$n) {
		nc <- b$nrows[i]*ncol(x)
		vv <- NULL
		if (expand) {
			v <- list(matrix(x[[1]]@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			v <- do.call(rbind, rep(v, halfway+1))
			for (k in 2:(1+halfway)) {
				v <- rbind(v,  matrix(x[[k]]@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			}
		} else if (pad) {
			v <- matrix(padvalue, ncol=b$nrows[i]*ncol(x), nrow=nread*halfway)
			for (k in 1:(1+halfway)) {
				v <- rbind(v,  matrix(x[[k]]@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=nc))
			}
		} else {
			v <- lapply(1:w[3],
				function(k) {
					y <- x[[k]]
					z <- y@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
					y <- messages(y, "focal3D")
					matrix(z, ncol=nc)
				})
				
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
					v <- rbind(v, matrix(x[[k]]@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt), ncol=ncol(v)))
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
				mv <- readValues(x[[j]], b$row[i], b$nrows[i])
				if (na.only) {
					k <- !is.na(mv)
				} else {
					k <- is.na(mv)
				}
				vout[k] <- mv[k]
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
	v <- x@cpp$focalValues(w, fillvalue, max(0, trunc(nrow(x)/2)), 1, opt)[1:prod(w)]
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
	b <- writeStart(out, filename, overwrite, n=msz+2, sources=sources(x), wopt=wopt)

	nc <- ncol(out)
	for (i in 1:b$n) {
		vv <- NULL
		for (j in 1:nl) {
			if (nl > 1) {
				v <- x[[j]]@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
			} else {
				v <- x@cpp$focalValues(w, fillvalue, b$row[i]-1, b$nrows[i], opt)
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


.getRegFun <- function(fun, weighted=FALSE, wopt, nmsx, nl, na.rm=FALSE, intercept=TRUE, ...) {

	ols <- function(x, y, ...) {
		if (any(is.na(x)) || any(is.na(y)) || (NROW(y) < (NCOL(x) + 1))) { 
			return(rep(NA, NCOL(x)+1))
		}
		stats::.lm.fit(cbind(1, x), y)$coefficients
	}

	ols_noi <- function(x, y, ...) {
		if (any(is.na(x)) || any(is.na(y)) || (NROW(y) < (NCOL(x)))) { 
			return(rep(NA, NCOL(x)))
		}
		stats::.lm.fit(as.matrix(x), y)$coefficients
	}

	ols_narm <- function(x, y, ...) {
		v <- stats::na.omit(cbind(y, x))
		if (nrow(v) < (NCOL(x) + 1)) {
			return( cbind(rep(NA, NCOL(x)+1)) )
		}
		stats::.lm.fit(cbind(1, v[,-1]), v[,1])$coefficients
	}
	
	ols_noi_narm <- function(x, y, ...) {
		v <- stats::na.omit(cbind(y, x))
		if (nrow(v) < (NCOL(x))) {
			return( cbind(rep(NA, NCOL(x))) )
		}
		stats::.lm.fit(v[,-1,drop=FALSE], v[,1])$coefficients
	}

	weighted_ols <- function(x, y, weights, ...) {
		if (any(is.na(x)) || any(is.na(y)) || any(is.na(weights)) || (NROW(y) < (NCOL(x) + 1))) { 
			return(rep(NA, NCOL(x)+1))
		}
		stats::lm.wfit(cbind(1, x), y, weights)$coefficients		
	}

	weighted_ols_noi <- function(x, y, weights, ...) {
		if (any(is.na(x)) || any(is.na(y)) || any(is.na(weights)) || (NROW(y) < (NCOL(x) + 1))) { 
			return(rep(NA, NCOL(x)))
		}
		stats::lm.wfit(as.matrix(x), y, weights)$coefficients		
	}

	weighted_ols_narm <- function(x, y, weights, ...) {
		v <- stats::na.omit(cbind(y, weights, x))
		if (nrow(v) < (NCOL(x) + 1)) {
			return(rep(NA, NCOL(x)+1))
		}
		stats::lm.wfit(cbind(1, v[,-c(1:2)]), v[,1], v[,2])$coefficients		
	}	

	weighted_ols_noi_narm <- function(x, y, weights, ...) {
		v <- stats::na.omit(cbind(y, weights, x))
		if (nrow(v) < (NCOL(x))) {
			return(rep(NA, NCOL(x)))
		}
		stats::lm.wfit(v[,-c(1:2), drop=FALSE], v[,1], v[,2])$coefficients		
	}	


	fun <- tolower(fun[1])
	if (fun != "ols") {
		return(list(fun=fun, wopt=wopt))
	}

	intercept <- isTRUE(intercept)
	
	if (intercept) {
		if (weighted) {
			if (na.rm) {
				fun <- weighted_ols_narm
			} else {
				fun <- weighted_ols				
			}
		} else {
			if (na.rm) {
				fun = ols_narm
			} else {
				fun = ols				
			}
		}
		if (is.null(wopt$names )) {
			wopt$names <- c("intercept", nmsx[-1])
		}
	} else {
		if (weighted) {
			if (na.rm) {
				fun <- weighted_ols_noi_narm
			} else {
				fun <- weighted_ols_noi				
			}
		} else {
			if (na.rm) {
				fun = ols_noi_narm
			} else {
				fun = ols_noi				
			}
		}
		if (is.null(wopt$names )) {
			wopt$names <- nmsx[-1]
		}
		nl = nl-1
	}

	list(fun=fun, wopt=wopt, nl=nl, intercept=intercept, na.rm=na.rm)
}

setMethod("focalReg", signature(x="SpatRaster"),
function(x, w=3, fun="ols", ..., fillvalue=NA, filename="", overwrite=FALSE, wopt=list()) {

	nl <- nlyr(x)
	if (nl < 2) error("focalReg", "x must have at least 2 layers")

	if (!is.numeric(w)) {
		error("focalReg", "w should be numeric vector or matrix")
	}

	weighted <- FALSE
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		m[m==0] <- NA
		test <- stats::na.omit(m)
		if (length(test) == 0) {
			error("focalReg", "all values in w are NA and/or zero")
		}
		if (any(test != 1)) {
			weighted <- TRUE
			message("the focal values are used as weights")
		} 
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}

	msz <- prod(w)
	if (msz < 2) {
		error("the effective weight matrix must have positive dimensions, and at least one must be > 1")
	}

	hasnam <- FALSE
	if (any(is.na(m))) {
		hasnam <- TRUE
		k <- !is.na(m)
		msz <- sum(k)
		weights <- m[k]
	} else if (weighted) {
		weights <- m	
	}

	if (is.character(fun)) {
		funopt <- .getRegFun(fun, weighted, wopt, names(x), nlyr(x), ...) 
		fun <- funopt$fun
		wopt <- funopt$wopt
		outnl <- funopt$nl
	} else {
		# need to test
		#outnl <- 
	}
	out <- rast(x, nlyr=outnl)
	
	b <- writeStart(out, filename, n=msz+2, sources=sources(x), wopt=wopt)
	ry <- x[[1]]
	rx <- x[[-1]]

	if (nl == 2) {
		for (i in 1:b$n) {
			Y <- focalValues(ry, w, b$row[i], b$nrows[i], fillvalue)
			X <- focalValues(rx, w, b$row[i], b$nrows[i], fillvalue)
			if (hasnam) {
				Y <- Y[,k,drop=FALSE]
				X <- X[,k,drop=FALSE]
			}
			if (weighted) {
				v <- t(sapply(1:nrow(Y), function(i) fun(X[i,], Y[i,], weights, ...)))			
			} else {
				v <- t(sapply(1:nrow(Y), function(i) fun(X[i,], Y[i,], ...)))
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	} else {

		for (i in 1:b$n) {
			Y <- focalValues(ry, w, b$row[i], b$nrows[i], fillvalue)
			if (hasnam) {
				Y <- Y[,k,drop=FALSE]
			}
			X <- list()
			for (j in 1:(nl-1)) {
				X[[j]] <- focalValues(rx[[j]], w, b$row[i], b$nrows[i], fillvalue)
				if (hasnam) {
					X[[j]] <- X[[j]][,k]
				}
			}
			v <- list()
			for (p in 1:nrow(Y)) {
				xlst <- list()
				for (j in 1:(nl-1)) {
					xlst[[j]] <- X[[j]][p,]
				}
				pX <- do.call(cbind, xlst)
				if (weighted) {
					v[[p]] <- fun(pX, Y[p,], weights)
				} else {
					v[[p]] <- fun(pX, Y[p,])
				}
			}
			v <- t(do.call(cbind, v))
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	}

	out <- writeStop(out)
	return(out)
}
)


setMethod("focalPairs", signature(x="SpatRaster"),
function(x, w=3, fun, ..., fillvalue=NA, filename="", overwrite=FALSE, wopt=list()) {

	pearson <- function(x, y, ...) { 
		.pearson(x, y, FALSE)
	}
	pearson_narm <- function(x, y, ...) { 
		.pearson(x, y, TRUE)
	}

	weighted_pearson <- function(x, y, weights, ...) { 
		.weighted_pearson(x, y, weights, FALSE)
	}
	weighted_pearson_narm <- function(x, y, weights, ...) { 
		.weighted_pearson(x, y, weights, TRUE)
	}
	
	nl <- nlyr(x)
	if (nl < 2) error("focalPairs", "x must have at least 2 layers")

	if (!is.numeric(w)) {
		error("focalPairs", "w should be numeric vector or matrix")
	}
	weighted <- FALSE
	if (is.matrix(w)) {
		m <- as.vector(t(w))
		m[m==0] <- NA
		test <- stats::na.omit(m)
		if (length(test) == 0) {
			error("focalPairs", "all values in w are NA and/or zero")
		}
		if (any(test != 1)) {
			weighted <- TRUE
			message("the focal values are used as weights")
		} 
		w <- dim(w)
	} else {
		w <- rep_len(w, 2)
		stopifnot(all(w > 0))
		m <- rep(1, prod(w))
	}
	msz <- prod(w)

	hasnam <- FALSE
	if (any(is.na(m))) {
		hasnam <- TRUE
		k <- !is.na(m)
		msz <- sum(k)
		weights <- m[k]
	} else if (weighted) {
		weights <- m	
	}

	if (msz < 2) {
		error("the effective weight matrix must have positive dimensions, and at least one must be > 1")
	}

	if (is.character(fun)) {
		narm <- isTRUE(list(...)$na.rm)
		fun <- tolower(fun[1])
		if (fun == "pearson") {
			if (weighted) {
				if (narm) {
					fun = weighted_pearson_narm
				} else {
					fun = weighted_pearson
				}
			} else {
				if (narm) {
					fun = pearson_narm
				} else {
					fun = pearson
				}
			}
		}
	}
	
	if (weighted) {
		test <- try(do.call(fun, list(1:prod(w), prod(w):1, weights=rep(1, prod(w)), ...)))
		if (inherits(test, "try-error")) {
			error("focalPairs", "'fun' does not work. Does it have a 'weights' argument?")
		}
	} else {
		test <- try(do.call(fun, list(1:prod(w), prod(w):1, ...)))
		if (inherits(test, "try-error")) {
			error("focalPairs", "'fun' does not work. Does it have two arguments (one for each layer)")
		}
	}
	if (is.null(wopt$names )) {
		wopt$names <- colnames(test)
	}

	outnl <- (nlyr(x) - 1) * length(test)
	out <- rast(x, nlyr=outnl)

	b <- writeStart(out, filename, n=msz+2, sources=sources(x), wopt=wopt)

	for (i in 1:b$n) {
		v <- list()
		Y <- focalValues(x[[1]], w, b$row[i], b$nrows[i], fillvalue)
		if (hasnam) {
			Y <- Y[,k,drop=FALSE]
		}
		for (j in 2:nlyr(x)) {
			X <- Y
			Y <- focalValues(x[[j]], w, b$row[i], b$nrows[i], fillvalue)
			if (hasnam) {
				Y <- Y[,k,drop=FALSE]
			} 
			if (weighted) {
				v[[j-1]] <- t(sapply(1:nrow(Y), function(i) fun(X[i,], Y[i,], weights=weights, ...)))
			} else {
				v[[j-1]] <- t(sapply(1:nrow(Y), function(i) fun(X[i,], Y[i,], ...)))
			}
		}
		v <- do.call(cbind, v)
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)




	# ..ols <- function(x, y, ...) {
		# v <- cbind(y, x)
		# if (any(is.na(v))) return( cbind(rep(NA, NCOL(x)+1)) )
		# X <- cbind(1, v[,-1])
		# XtX <- t(X) %*% X
		# if (det(XtX) == 0) {
			# return(rep(NA, NCOL(y)+1))
		# }
		# invXtX <- solve(XtX) %*% t(X)
		# invXtX %*% v[,1]
	# }

	# ..ols_noi <- function(x, y, ...) {
		# v <- cbind(y, x)
		# if (any(is.na(v))) return( cbind (rep(NA, NCOL(x))) )
		# X <- v[,-1,drop=FALSE]
		# XtX <- t(X) %*% X
		# if (det(XtX) == 0) {
			# return(rep(NA, ncol(y)+1))
		# }
		# invXtX <- solve(XtX) %*% t(X)
		# invXtX %*% v[,1]
	# }

	# ..ols_narm <- function(x, y, ...) {
		# v <- na.omit(cbind(y, x))
		# if (nrow(v) < (NCOL(x) + 1)) {
			# return( cbind(rep(NA, NCOL(x)+1)) )
		# }
		# X <- cbind(1, v[,-1])
		# XtX <- t(X) %*% X
		# if (det(XtX) == 0) {
			# return(NA)
		# }
		# invXtX <- solve(XtX) %*% t(X)
		# invXtX %*% v[,1]
	# }
	
	# ..ols_noi_narm <- function(x, y, ...) {
		# v <- na.omit(cbind(y, x))
		# if (nrow(v) < NCOL(x)) {
			# return( cbind(rep(NA, NCOL(y))) )
		# }
		# X <- v[,-1,drop=FALSE]
		# XtX <- t(X) %*% X
		# if (det(XtX) == 0) {
			# return(NA)
		# }
		# invXtX <- solve(XtX) %*% t(X)
		# invXtX %*% v[,1]
	# }

	# ..weighted_ols <- function(x, y, weights, ...) {
		# if (any(is.na(x)) || any(is.na(y))) { 
			# return(rep(NA, NCOL(x)+1))
		# }
		# stats::coefficients(stats::glm(y~x, weights=weights))
	# }

	# ..weighted_ols_noi <- function(x, y, weights, ...) {
		# if (any(is.na(x)) || any(is.na(y))) { 
			# return(rep(NA, NCOL(x)))
		# }	
		# stats::coefficients(stats::glm(y ~ -1 + ., weights=weights))
	# }

	# ..weighted_ols_narm <- function(x, y, weights, ...) {
		# v <- na.omit(data.frame(y=y, x, weights=weights))
		# if (nrow(v) < (NCOL(x) + 1)) {
			# return(rep(NA, NCOL(x)+1))
		# }
		# weights <- v$weights
		# v$weights <- NULL
		# stats::coefficients(stats::glm(y ~ ., data=v, weights=weights))
	# }	

	# ..weighted_ols_noi_narm <- function(x, y, weights, ...) {
		# v <- na.omit(data.frame(y=y, x, weights=weights))
		# if (nrow(v) < (NCOL(x))) {
			# return(rep(NA, NCOL(x)))
		# }
		# weights <- v$weights
		# v$weights <- NULL
		# stats::coefficients(stats::glm(y ~ -1 + ., data=v, weights=weights))
	# }	
