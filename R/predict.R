# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3


.runModel <- function(model, fun, d, nl, const, na.rm, index, ...) {
	if (! is.null(const)) {
		d <- cbind(d, const[1])
	} 
	if (na.rm) {
		n <- nrow(d)
		i <- rowSums(is.na(d)) == 0
		d <- d[i,,drop=FALSE]
		if (ncol(d) > 0) {
			r <- fun(model, d, ...)
			if (is.factor(r)) {
				r <- as.integer(r)
			} else if (is.data.frame(r)) {
				r <- sapply(r, as.numeric)
			}
			r <- as.matrix(r)
			if (!all(i)) {
				m <- matrix(NA, nrow=nl*n, ncol=ncol(r))
				m[i,] <- r
				colnames(m) <- colnames(r)
				r <- m
			}
		} else {
			r <- matrix(NA, nrow=nl*n, ncol=1)
		}
	} else {
		r <- fun(model, d, ...)
		if (is.factor(r)) {
			r <- as.integer(r)
		} else if (is.data.frame(r)) {
			r <- sapply(r, as.numeric)
		}
		r <- as.matrix(r)
	}
	if (!is.null(index)) {
		r <- r[, index,drop=FALSE]
	}
	r
}



.getFactors <- function(m, facts=NULL, lyrnms) {

	if (!is.null(facts)) {
		stopifnot(is.list(factors))
		f <- names(factors)
		if (any(trimws(f) == "")) {
			stop("all factors must be named")
		}
	} else if (inherits(m, "randomForest")) {
		f <- names(which(sapply(m$forest$xlevels, max) != "0"))
		if (length(f) > 0) { 
			factors <- m$forest$xlevels[f]
		}
	} else if (inherits(m, "gbm")) {
		dafr <- m$gbm.call$dataframe 
		i <- sapply(dafr, is.factor)
		if (any(i)) {
			j <- which(i)
			factors <- list()
			for (i in 1:length(j)) {
				factors[[i]] <- levels(dafr[[ j[i] ]])
			}
			names(factors) <- colnames(dafr)[j]
		}
	} else { #glm and others
		try(factors <- m$xlevels, silent=TRUE)
	}
	if (!all(names(factors) %in% lyrnms)) {
		ff <- f[!(f %in% lyrnms)]
		stop(paste("factor name(s):", paste(ff, collapse=", "), " not in layer names"))
	}
	factors
}
	
setMethod("predict", signature(object="SpatRaster"), 
	function(object, model, fun=predict, ..., factors=NULL, const=NULL, na.rm=FALSE, index=NULL, filename="", overwrite=FALSE, wopt=list()) {

		nms <- names(object)
		if (length(unique(nms)) != length(nms)) {
			tab <- table(nms)
			stop('duplicate names in SpatRaster: ', tab[tab>1])
		}
		
		#factors should come with the SpatRaster
		#haveFactor <- FALSE
		#if (!is.null(factors)) {
		#	factors <- .getFactors(model, factors, nms)
		#	fnames <- names(f)
		#	haveFactor <- TRUE
		#}
		
		nl <- 1
		nc <- ncol(object)
		tomat <- FALSE
		d <- readValues(object, round(0.5*nrow(object)), 1, 1, min(nc,500), TRUE, TRUE)

		r <- .runModel(model, fun, d, nl, const, na.rm, index, ...)
		nl <- ncol(r)		
		out <- rast(object, nlyr=nl)
		cn <- colnames(r)
		if (length(cn) == nl) names(out) <- make.names(cn, TRUE)
		
		readStart(object)
		b <- writeStart(out, filename, overwrite, wopt)
		for (i in 1:b$n) {
			d <- readValues(object, b$row[i], b$nrows[i], 1, nc, TRUE, TRUE)
			r <- .runModel(model, fun, d, nl, const, na.rm, index, ...)
			writeValues(out, r, b$row[i], b$nrows[i])
		}
		readStop(object)
		out <- writeStop(out)
		return(out)
	}
)


