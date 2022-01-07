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
			depf <- deparse(fun)
			test1 <- isTRUE(try( depf[2] == 'UseMethod(\"mean\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "mean", silent=TRUE))
			if (test1 | test2) { 
				fun <- "mean" 
			}
			test1 <- isTRUE(try( depf[2] == 'UseMethod(\"median\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "median", silent=TRUE))
			if (test1 | test2) { 
				fun <- "median" 
			}
			test1 <- isTRUE(try( depf[1] == "function (x, na.rm = FALSE) ", silent=TRUE))
			test2 <- isTRUE(try( depf[2] == "sqrt(var(if (is.vector(x) || is.factor(x)) x else as.double(x), ", silent=TRUE))
			test3 <- isTRUE(try( depf[3] == "    na.rm = na.rm))", silent=TRUE))
			if (test1 && test2 && test3) { 
				fun <- "sd"
			}
		} 
	}
	return(fun)
}


setMethod("aggregate", signature(x="SpatRaster"), 
function(x, fact=2, fun="mean", ..., cores=1, filename="", overwrite=FALSE, wopt=list())  {

	fun <- .makeTextFun(fun)
	toc <- FALSE
	if (class(fun) == "character") { 
		if (fun %in% c("sum", "mean", "min", "max", "median", "modal", "sd", "sdpop")) {
			toc <- TRUE
		} else {
			fun <- match.fun(fun) 
		}
	} else {
		fun <- match.fun(fun) 
	}
	if (!hasValues(x)) { toc = TRUE }
	if (toc) {
		#	fun="mean", expand=TRUE, na.rm=TRUE, filename=""
		narm <- isTRUE(list(...)$na.rm)
		opt <- spatOptions(filename, overwrite, wopt=wopt)
		x@ptr <- x@ptr$aggregate(fact, fun, narm, opt)
		return (messages(x, "aggregate"))
	} else {
		out <- rast(x)
		nl <- nlyr(out)
		opt <- spatOptions()
		out@ptr <- out@ptr$aggregate(fact, "sum", TRUE, opt)
		out <- messages(out, "aggregate")

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
		on.exit(readStop(x))
		ignore <- writeStart(out, filename, overwrite, wopt=wopt)
		if (doPar) {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc)
				v <- x@ptr$get_aggregates(v, b$nrows[i], dims)
				v <- parallel::parSapply(cls, v, fun, ...)
				if (length(v) != outnr[i] * prod(dims[5:6])) {
					error("aggregate", "this function does not return the correct number of values")
				}
				writeValues(out, v, outrows[i], outnr[i])
			}
		} else {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc)
				v <- x@ptr$get_aggregates(v, b$nrows[i], dims)
				v <- sapply(v, fun, ...)
				if (length(v) != outnr[i] * prod(dims[5:6])) {
					error("aggregate", "this function does not return the correct number of values")
				}
				writeValues(out, v, outrows[i], outnr[i])
			}
		}
		out <- writeStop(out)
		messages(out, "aggregate")
	}
}
)


.agg_uf <- function(i) {
	u <- unique(i)
	if (length(u) == 1) { u } else { NA	}
}


aggregate_attributes <- function(d, by, fun=NULL, ...) {
	i <- sapply(d, is.numeric)
	i[colnames(d) %in% by] <- FALSE
	j <- 1:length(by)
	da <- db <- NULL
	if (!is.null(fun)) {
		if (any(i)) {
			if (is.character(fun)) {
				f <- match.fun(fun)
				da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], f)
				names(da)[-j] <- paste0(fun, "_", names(da)[-j])
			} else {
				da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], fun)
				names(da)[-j] <- paste0("agg_", names(da)[-j])
			}
		}
	} else {
		i[] <- FALSE
	}
	i[colnames(d) %in% by] <- TRUE
	if (any(!i)) {
		db <- aggregate(d[, !i,drop=FALSE], d[, by, drop=FALSE], .agg_uf)
		db <- db[, colSums(is.na(db)) < nrow(db), drop=FALSE]
		if (NCOL(da)>1) {
			da <- merge(da, db, by=by)
		} else {
			da <- db
		}
	}

	dn <- aggregate(d[, by,drop=FALSE], d[, by, drop=FALSE], length)
	colnames(dn)[2] = "agg_n"
	if (NCOL(da)>1) {
		dn <- merge(da, dn, by=by)
	}
	dn
}


setMethod("aggregate", signature(x="SpatVector"),
	function(x, by=NULL, dissolve=TRUE, fun="mean", ...) {
		if (length(by) > 1) {
			# to be fixed
			error("aggregate", "this method can only aggregate by one variable")
		}
		if (is.numeric(by[1])) {
			 i <- round(by)
			 if ((i > 0) & (i <= ncol(x))) {
				 by <- names(x)[i]
			 } else {
				 error("aggregate", "invalid column number supplied: ", by)
			 }
		}

		if (is.null(by)) {
			x$aggregate_by_variable = 1;
			x@ptr <- x@ptr$aggregate("aggregate_by_variable", dissolve)
			x$aggregate_by_variable = NULL;
		} else {
			d <- as.data.frame(x)
			x@ptr <- x@ptr$aggregate(by, dissolve)
			a <- aggregate_attributes(d, by, fun)
			if (any(is.na(d[[by]]))) {
				# because NaN and NA are dropped
				i <- nrow(a)+(1:2)
				a[i,] <- c(NA, NaN)
				i <- match(a[[by]], x[[by,drop=TRUE]])
				i <- i[!is.na(i)]
			} else {
				i <- match(a[[by]], x[[by,drop=TRUE]])			
			}
			values(x) <- a[i,]
		}
		x
	}
)

# setMethod("aggregate", signature(x="SpatVector"),
	# function(x, by=NULL, dissolve=TRUE, fun="mean", ...) {
		# gt <- geomtype(x)
		# if (length(by) > 1) {
			# error("aggregate", "this method can only aggregate by one variable")
		# }
		# x <- methods::as(x, "Spatial")
		# if (is.numeric(by[1])) {
			# i <- round(by)
			# if ((i > 0) & (i <= ncol(x))) {
				# by <- names(x)[i]
			# } else {
				# error("aggregate", "invalid column number supplied: ", by)
			# }
		# }
		# r <- aggregate(x, by=by, dissolve=dissolve, ...)
		# if (!missing(fun) && !missing(by)) {
			# if (.hasSlot(x, "data")) {
				# d <- x@data
				# i <- sapply(d, is.numeric)
				# i[colnames(d) %in% by] <- FALSE
				# j <- 1:length(by)
				# if (any(i)) {
					# if (is.character(fun)) {
						# f <- match.fun(fun)
						# da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], f)
						# names(da)[-j] <- paste0(fun, "_", names(da)[-j])
					# } else {
						# da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], fun)
						# names(da)[-j] <- paste0("agg_", names(da)[-j])
					# }
					# r <- merge(r, da, by)
				# }
				# i[colnames(d) %in% by] <- TRUE
				# if (any(!i)) {
					# db <- aggregate(d[, !i,drop=FALSE], d[, by, drop=FALSE], .agg_uf)
					# db <- db[, colSums(is.na(db)) < nrow(db), drop=FALSE]
					# if (ncol(db) > 1) {
						# r <- merge(r, db, by)
					# }
				# }
				# dn <- aggregate(d[, by,drop=FALSE], d[, by, drop=FALSE], length)
				# colnames(dn)[2] = "agg_n"
				# r <- merge(r, dn, by)
			# }
		# }
		# vect(r)
	# }
# )

