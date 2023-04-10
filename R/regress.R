# Author: Robert J. Hijmans
# Date : April 2023
# Version 1.0
# License GPL v3



.get_conY_MM <- function(formula, Y, na.rm=FALSE, nl) {
	formula <- eval(as.formula(formula))
	nas <- eval(rep(NA, nl))
	if (na.rm) {
		ols <- function(x, ...) {
			d <- na.omit(data.frame(y=Y, x=x))
			mm <- model.matrix(formula, data=d)
			if (nrow(mm < (nl))) {
				return(nas)
			}
			stats::.lm.fit(mm, d$y)$coefficients
		}
	} else {
		ols <- function(x, ...) {
			if (any(is.na(x))) {
				return(nas)
			}
			d <- data.frame(y=Y, x=x)
			mm <- model.matrix(formula, data=d)
			stats::.lm.fit(mm, Y)$coefficients
		}
	}
	ols
}


setMethod("regress", signature(x="SpatRaster", y="numeric"),
function(x, y, formula=y~x, ..., cores=1, filename="", overwrite=FALSE, wopt=list()) {

	if (any(is.na(y))) {
		error("regress", "y cannot have NAs")
	}
	formula <- as.formula(formula)
	dat <- data.frame(y=y, x=1)
	mm <- model.matrix(formula, data=dat)
	outnl <- ncol(mm)		
	regfun <- .get_conY_MM(formula, Y=y, nl=outnl, ...)
	
	out <- rast(x)
	nlyr(out) <- outnl
	names(out) <- colnames(mm)
	nc <- ncol(x)
	readStart(x)
	on.exit(readStop(x))
	nl <- nlyr(x)

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
			r <- parallel::parRapply(cores, v, regfun, chunk.size=icsz)
			if (nlyr(out) > 1) {
				r <- matrix(r, ncol=nlyr(out), byrow=TRUE)
			}
			writeValues(out, r, b$row[i], b$nrows[i])
		}
	} else {
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
			r <- apply(v, 1, regfun)
			writeValues(out, t(r), b$row[i], b$nrows[i])
		}
	}
	writeStop(out)
})


setMethod("regress", signature(x="SpatRaster", y="missing"),
function(x, y, formula=y~x, ..., filename="", overwrite=FALSE, wopt=list()) {
	tm <- time(x)
	if (any(is.na(tm))) {
		y <- 1:nlyr(x)
	} else {
		y <- as.numeric(tm)
	}
	regress(x, y=y, formula=formula, ..., filename=filename, overwrite=overwrite, wopt=wopt)
})


setMethod("regress", signature(x="SpatRaster", y="SpatRaster"),
function(x, y, formula=y~x, ..., filename="", overwrite=FALSE, wopt=list()) {


})


