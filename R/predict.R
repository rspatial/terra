# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3

parfun <- function(cls, d, fun, model, ...) {
	nr <- nrow(d)
	nc <- length(cls)
	s <- split(d, rep(1:nc, each=ceiling(nr/nc), length.out=nr))
	p <- parallel::clusterApply(cls, s, function(i, ...) fun(model, i, ...), ...)
	if (!is.null(dim(p[[1]]))) {
		do.call(rbind, p)
	} else {
		unlist(p)
	}
}


.runModel <- function(model, fun, d, nl, const, na.rm, index, cores, ...) {
	doPar <- !is.null(cores)
	if (!is.data.frame(d)) {
		d <- data.frame(d)
	}
	if (!is.null(const)) {
		for (i in 1:ncol(const)) {
			d <- cbind(d, const[,i,drop=FALSE])
		}
	}
	if (na.rm) {
		n <- nrow(d)
		i <- rowSums(is.na(d)) == 0
		d <- d[i,,drop=FALSE]
		if (nrow(d) > 0) {
			if (doPar) {
				r <- parfun(cores, d, fun, model, ...)
			} else {
				r <- fun(model, d, ...)
			}
			if (is.list(r)) {
				r <- as.data.frame(r)
				# data.frame(lapply) instead of sapply to catch a one-row case
				r <- data.frame(lapply(r, as.numeric))
			} else if (is.factor(r)) {
				r <- as.integer(r)
			} 
			r <- as.matrix(r)
			if (!all(i)) {
				m <- matrix(NA, nrow=n, ncol=ncol(r))
				m[i,] <- r
				colnames(m) <- colnames(r)
				r <- m
			}
		} else {
			if (!is.null(index)) {
				r <- matrix(NA, nrow=nl*n, ncol=max(index))
			} else {
				r <- matrix(NA, nrow=nl*n, ncol=1)
			}
		}
	} else {
		if (doPar) {
			r <- parfun(cores, d, fun, model, ...)
		} else {
			r <- fun(model, d, ...)
		}
		if (is.list(r)) {
			r <- as.data.frame(lapply(r, as.numeric))
		} else if (is.factor(r)) {
			r <- as.integer(r)
		} else if (is.data.frame(r)) {
			r <- sapply(r, as.numeric)
		}
		r <- as.matrix(r)
	}
	if (inherits(model, "gstat")) {
		if (ncol(r) > 2) {
			nr <- max(nrow(d), 5)
			xy <- as.matrix(d[1:nr,1:2])
			if (all(xy == r[1:nr, 1:2])) {
				r <- r[,-c(1:2)]   # x, y
			}
		}
	}
	if (!is.null(index)) {
		r <- r[, index, drop=FALSE]
	}
	r
}


.getFactors <- function(model, fun, d, nl, const, na.rm, index, ...) {
	if (!is.data.frame(d)) {
		d <- data.frame(d)
	}
	if (! is.null(const)) {
		for (i in 1:ncol(const)) {
			d <- cbind(d, const[,i,drop=FALSE])
		}
	}
	if (na.rm) {
		n <- nrow(d)
		i <- rowSums(is.na(d)) == 0
		d <- d[i,,drop=FALSE]
	}
	if (nrow(d) > 0) {
		r <- fun(model, d, ...)
	} else {
		return(NULL)
	}

	if (inherits(model, "gstat")) {
		nr <- max(nrow(d), 5)
		xy <- d[1:nr,1:2]
		if (all(xy == r[1:nr, 1:2])) {
			r <- r[,-c(1:2)]   # x, y
		}
	}

	if (is.factor(r)) {
		levs <- levels(r)
		data.frame(value=1:length(levs), class=levs)
	} else if (is.list(r) || is.data.frame(r)) {
		r <- as.data.frame(r)
		out <- sapply(r, levels)
		for (i in 1:length(out)) {
			if (!is.null(out[[i]])) {
				out[[i]] <- data.frame(value=1:length(out[[i]]), label=out[[i]])
			}
		}
		out
	} else {
		NULL
	}
}


find_dims <- function(object, model, nc, fun, const, na.rm, index, ...) {
	nr <- nrow(object)
	nl <- 1
	testrow <- round(0.51*nr)
	rnr <- 1
	if (nc==1) rnr <- min(nr, 20) - testrow + 1
	d <- readValues(object, testrow, rnr, 1, nc, TRUE, TRUE)
	cn <- NULL
	levs <- NULL
	if (!is.null(index)) {
		nl <- length(index)
		r <- .runModel(model, fun, d, nl, const, na.rm, index, cores=NULL, ...)
		rdim <- dim(r)
		if (is.null(rdim)) {
			cn <- ""
		} else {
			cn <- colnames(r)
		}
	} else {
		allna <- FALSE
		if (na.rm) {
			allna <- all(nrow(stats::na.omit(d)) == 0)
			if (allna) {
				testrow <- ceiling(testrow - 0.25*nr)
				d <- readValues(object, testrow, rnr, 1, nc, TRUE, TRUE)
				allna <- all(nrow(stats::na.omit(d)) == 0)
			}
			if (allna) {
				testrow <- floor(testrow + 0.5*nr)
				if ((testrow + rnr) > nr) rnr = nr - testrow + 1
				d <- readValues(object, testrow, rnr, 1, nc, TRUE, TRUE)
				allna <- all(nrow(stats::na.omit(d)) == 0)
			}
			if (allna && (ncell(object) < 1000)) {
				d <- readValues(object, 1, nr, 1, nc, TRUE, TRUE)
				allna <- all(nrow(stats::na.omit(d)) == 0)
				#if (allna) {
				#	error("predict", "all predictor values are NA")
				#}
			}
			if (allna) {
				d <- spatSample(object, min(1000, ncell(object)), "regular", warn=FALSE)
				allna <- all(nrow(stats::na.omit(d)) == 0)
			}
			if (allna) {
				d[] <- stats::runif(prod(dim(d)))
			}
		}
		r <- .runModel(model, fun, d, nl, const, na.rm, index, cores=NULL, ...)
		if (ncell(object) > 1) {
			rdim <- dim(r)
			if (is.null(rdim)) {
				nl <- 1
				cn <- ""
			} else {
				if (isTRUE(any(rdim == 1))) {
					nl <- 1
					cn <- colnames(r)[1]
				} else {
					nl <- ncol(r)
					cn <- colnames(r)
				}
			}
		} else {
			nl <- length(r)
		}
		levs <- .getFactors(model, fun, d, nl, const, na.rm, index, ...)
	}
	out <- rast(object, nlyrs=nl)
	if (!all(sapply(levs, is.null))) levels(out) <- levs
	if (length(cn) == nl) names(out) <- make.names(cn, TRUE)
	out
}


setMethod("predict", signature(object="SpatRaster"),
	function(object, model, fun=predict, ..., const=NULL, na.rm=FALSE, index=NULL, cores=1, cpkgs=NULL, filename="", overwrite=FALSE, wopt=list()) {

		nms <- names(object)
		if (length(unique(nms)) != length(nms)) {
			tab <- table(nms)
			error("predict", "duplicate names in SpatRaster: ", tab[tab>1])
		}

		nc <- ncol(object)
		#tomat <- FALSE
		readStart(object)
		on.exit(readStop(object))

		out <- find_dims(object, model, nc, fun, const, na.rm, index, ...)
		nl <- nlyr(out)
		
		doclust <- FALSE
		if (inherits(cores, "cluster")) {
			doclust <- TRUE
		} else if (cores > 1) {
			doclust <- TRUE
			cores <- parallel::makeCluster(cores)
			on.exit(parallel::stopCluster(cores), add=TRUE)
		}
		if (doclust) {
			parallel::clusterExport(cores, c("model", "fun"), environment())
			if (!is.null(cpkgs)) {
				parallel::clusterExport(cores, "cpkgs", environment())
				parallel::clusterCall(cores, function() for (i in 1:length(cpkgs)) {library(cpkgs[i], character.only=TRUE) })
			}
			export_args(cores, ..., caller="predict")
			
		} else {
			cores <- NULL
		}
		b <- writeStart(out, filename, overwrite, wopt=wopt, n=max(nlyr(out), nlyr(object))*4, sources=sources(object))
		for (i in 1:b$n) {
			d <- readValues(object, b$row[i], b$nrows[i], 1, nc, TRUE, TRUE)
			r <- .runModel(model, fun, d, nl, const, na.rm, index, cores=cores, ...)
			if (prod(NROW(r), NCOL(r)) != prod(b$nrows[i], nc, nl)) {
				msg <- "the number of values returned by 'fun' (model predict function) does not match the input."
				if (!na.rm) msg <- paste(msg, "Try na.rm=TRUE?")
				error("predict", msg)
			}
			writeValues(out, r, b$row[i], b$nrows[i])
		}
		writeStop(out)
#		return(out)
	}
)

