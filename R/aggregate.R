# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


.makeTextFun <- function(fun) {
	if (!inherits(fun, "character")) {
		fun <- match.fun(fun)
		if (is.primitive(fun)) {
			test <- try(deparse(fun)[[1]], silent=TRUE)
			if (test == '.Primitive(\"sum\")') { fun <- 'sum' 
			} else if (test == '.Primitive(\"min\")') { fun <- 'min' 
			} else if (test == '.Primitive(\"max\")') { fun <- 'max' 
			}
		} else {
			test1 <- isTRUE(try( deparse(fun)[2] == 'UseMethod(\"mean\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "mean", silent=TRUE))
			if (test1 | test2) { 
				fun <- "mean" 
			}
			test1 <- isTRUE(try( deparse(fun)[2] == 'UseMethod(\"median\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "median", silent=TRUE))
			if (test1 | test2) { 
				fun <- "median" 
			}
		} 
	}
	return(fun)
}


setMethod("aggregate", signature(x="SpatRaster"), 
function(x, fact=2, fun="mean", ..., cores=1, filename="", overwrite=FALSE, wopt=list())  {

	fun <- .makeTextFun(match.fun(fun))
	toc <- FALSE
	if (class(fun) == "character") { 
		if (fun %in% c("sum", "mean", "min", "max", "median", "modal")) {
			toc <- TRUE
		}
	}
	if (!hasValues(x)) { toc = TRUE }
	if (toc) {	
		#	fun="mean", expand=TRUE, na.rm=TRUE, filename=""
		narm <- isTRUE(list(...)$na.rm)
		opt <- .runOptions(filename, overwrite, wopt)	
		x@ptr <- x@ptr$aggregate(fact, fun, narm, opt)
		return (show_messages(x, "aggregate"))
	} else {
		out <- rast(x)
		nl <- nlyr(out)
		opt <- .getOptions()	
		out@ptr <- out@ptr$aggregate(fact, "sum", TRUE, opt)
		out <- show_messages(out, "aggregate")
		
		dims <- x@ptr$get_aggregate_dims(fact)
		b <- x@ptr$getBlockSize(4, opt$memfrac)		
		
		nr <- max(1, floor(b$nrows[1] / fact[1])) * fact[1]
		nrs <- rep(nr, floor(nrow(x)/nr))
		d <- nrow(x) - sum(nrs) 
		if (d > 0) nrs <- c(nrs, d)
		b$row <- c(0, cumsum(nrs))[1:length(nrs)] + 1
		b$nrows <- nrs
		b$n <- length(nrs)
		outnr <- ceiling(b$nrows / fact[1])
		outrows  <- c(0, cumsum(outnr))[1:length(outnr)] + 1	
		nc <- ncol(x)	
		
		if (cores > 1) {
			doPar <- TRUE
			cls <- parallel::makeCluster(cores)
			on.exit(parallel::stopCluster(cls))
			#f <- function(v, ...) parallel::parSapply(cls, v, fun, ...)
		} else {
			doPar <- FALSE
			#f <- function(v, ...) sapply(v, fun, ...)
		}

		readStart(x)
		ignore <- writeStart(out, filename, overwrite, wopt)
		if (doPar) {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc)
				v <- x@ptr$get_aggregates(v, b$nrows[i], dims)
				v <- parallel::parSapply(cls, v, fun, ...)
				if (length(v) != outnr[i] * prod(dims[5:6])) {
					stop("this function does not return the correct number of values")
				}
				writeValues(out, v, outrows[i], outnr[i])
			}	
		} else {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc)
				v <- x@ptr$get_aggregates(v, b$nrows[i], dims)
				v <- sapply(v, fun, ...)	
				if (length(v) != outnr[i] * prod(dims[5:6])) {
					stop("this function does not return the correct number of values")
				}
				writeValues(out, v, outrows[i], outnr[i])
			}	
		}
		readStop(x)
		out <- writeStop(out)
		show_messages(out, "aggregate")
	}
}
)


.agg_uf <- function(i) {
	u <- unique(i)
	if (length(u) == 1) { u } else { NA	}
}


## cheating --- using raster/sp/rgeos for now
## but improved over raster
setMethod("aggregate", signature(x="SpatVector"),
	function(x, by=NULL, dissolve=TRUE, fun="mean", ...) {
		#gt <- geomtype(x)
		if (length(by) > 1) {
			stop("this method can only aggregate by one variable")
		}
		x <- methods::as(x, "Spatial")
		if (is.numeric(by[1])) {
			i <- round(by)
			if ((i > 0) & (i <= ncol(x))) {
				by <- names(x)[i]
			} else {
				stop(paste("invalud column number supplied:", by))
			}
		}
		r <- aggregate(x, by=by, dissolve=dissolve, ...)
		if (!missing(fun) && !missing(by)) {
			if (.hasSlot(x, "data")) {
				d <- x@data
				i <- sapply(d, is.numeric)
				i[colnames(d) %in% by] <- FALSE
				j <- 1:length(by)
				if (any(i)) {
					if (is.character(fun)) {
						f <- match.fun(fun)
						da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], f)				
						names(da)[-j] <- paste0(fun, "_", names(da)[-j])
					} else {
						da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], fun)
						names(da)[-j] <- paste0("agg_", names(da)[-j])
					}
					r <- merge(r, da, by)
				}
				i[colnames(d) %in% by] <- TRUE
				if (any(!i)) {
					db <- aggregate(d[, !i,drop=FALSE], d[, by, drop=FALSE], .agg_uf)	
					db <- db[, colSums(is.na(db)) < nrow(db), drop=FALSE]
					if (ncol(db) > 1) {
						r <- merge(r, db, by)
					}
				}
				dn <- aggregate(d[, by,drop=FALSE], d[, by, drop=FALSE], length)
				colnames(dn)[2] = "agg_n"
				r <- merge(r, dn, by)
			}
		}
		vect(r)
	}
)

