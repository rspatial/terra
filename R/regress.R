# Author: Robert J. Hijmans
# Date : April 2023
# Version 1.0
# License GPL v3



.reg_constX <- function(formula, X, na.rm=FALSE, nl) {
	formula <- eval(as.formula(formula))
	nas <- eval(rep(NA, nl))
	d <- data.frame(x=X, y=X)
	mm <- eval(model.matrix(formula, data=d))
	if (na.rm) {
		ols <- function(y, ...) {
			m <- na.omit(cbind(y, mm))
			if (nrow(m) == 0) {
				return(nas)
			}
			stats::.lm.fit(m[,-1,drop=FALSE], m[,1])$coefficients
		}
	} else {
		ols <- function(y, ...) {
			if (any(is.na(y))) {
				return(nas)
			}
			stats::.lm.fit(mm, y)$coefficients
		}
	}
	ols
}



setMethod("regress", signature(y="SpatRaster", x="numeric"),
function(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...) {

	if (any(is.na(x))) {
		error("regress", "y cannot have NAs")
	}
	formula <- as.formula(formula)
	dat <- data.frame(x=x, y=1)
	mm <- model.matrix(formula, data=dat)
	outnl <- ncol(mm)		
	regfun <- .reg_constX(formula, X=x, na.rm=na.rm, nl=outnl)
	
	out <- rast(y)
	nlyr(out) <- outnl
	names(out) <- colnames(mm)
	nc <- ncol(y)
	readStart(y)
	on.exit(readStop(y))
	nl <- nlyr(y)

	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		doclust <- TRUE
	} else if (cores > 1) {
		doclust <- TRUE
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores), add=TRUE)
	}

	ncops <- nlyr(y) / nlyr(out)
	ncops <- ifelse(ncops > 1, ceiling(ncops), 1) * 4
	b <- writeStart(out, filename, overwrite, n=ncops, sources=sources(y), ...)

	if (doclust) {
		ncores <- length(cores)
		#export_args(cores, ...)
		for (i in 1:b$n) {
			v <- readValues(y, b$row[i], b$nrows[i], 1, nc, TRUE)
			icsz <- max(min(100, ceiling(b$nrows[i] / ncores)), b$nrows[i])
			r <- parallel::parRapply(cores, v, regfun, chunk.size=icsz)
			if (nlyr(out) > 1) {
				r <- matrix(r, ncol=nlyr(out), byrow=TRUE)
			}
			writeValues(out, r, b$row[i], b$nrows[i])
		}
	} else {
		for (i in 1:b$n) {
			v <- readValues(y, b$row[i], b$nrows[i], 1, nc, TRUE)
			r <- apply(v, 1, regfun)
			writeValues(out, t(r), b$row[i], b$nrows[i])
		}
	}
	writeStop(out)
})


setMethod("regress", signature(x="SpatRaster", y="missing"),
function(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...) {
	tm <- time(y)
	if (any(is.na(tm))) {
		x <- 1:nlyr(y)
	} else {
		x <- as.numeric(tm)
	}
	regress(y, x=x, formula=formula, na.rm=na.rm, cores=cores, filename=filename, overwrite=overwrite, ...)
})


setMethod("regress", signature(x="SpatRaster", y="SpatRaster"),
function(y, x, formula=y~x, ..., filename="", overwrite=FALSE, wopt=list()) {
	error("regress", "to be implemented")
})


