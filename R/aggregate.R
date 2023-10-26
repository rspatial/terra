# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3


.makeTextFun <- function(fun) {
	if (!inherits(fun, "character")) {
		fun <- match.fun(fun)
		if (is.primitive(fun)) {
			test <- try(deparse(fun)[[1]], silent=TRUE)
			test <- gsub('.Primitive\\(\"', "", test)
			test <- gsub('\")', "", test)
			if (test %in% c("sum", "min", "max", "prod", "any", "all")) return(test);
		} else {
			depf <- deparse(fun)
			test1 <- isTRUE(try( depf[2] == 'UseMethod(\"mean\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "mean", silent=TRUE))
			if (test1 | test2) return("mean")
			test1 <- isTRUE(try( depf[2] == 'UseMethod(\"median\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "median", silent=TRUE))
			if (test1 | test2) return("median")
			test1 <- isTRUE(try( depf[1] == "function (x, na.rm = FALSE) ", silent=TRUE))
			test2 <- isTRUE(try( depf[2] == "sqrt(var(if (is.vector(x) || is.factor(x)) x else as.double(x), ", silent=TRUE))
			test3 <- isTRUE(try( depf[3] == "    na.rm = na.rm))", silent=TRUE))
			if (test1 && test2 && test3) return("sd")
			if (isTRUE(try( fun@generic == "which.min", silent=TRUE))) return("which.min")
			if (isTRUE(try( fun@generic == "which.max", silent=TRUE))) return("which.max")
			if (isTRUE(all(depf[1] == deparse(base::which)[1]))) return("which")
			if (isTRUE(all(depf[1] == deparse(base::table)[1]))) return("table")
		}
	}
	return(fun)
}


setMethod("aggregate", signature(x="SpatRaster"),
function(x, fact=2, fun="mean", ..., cores=1, filename="", overwrite=FALSE, wopt=list())  {

	fun <- .makeTextFun(fun)
	toc <- FALSE
	if (inherits(fun, "character")) {
		if (fun %in% c("sum", "mean", "min", "max", "median", "modal","prod", "which.min", "which.max",
				"any", "all", "sd", "std", "sdpop")) {
			fun[fun == "sdpop"] <- "std"
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
		x@cpp <- x@cpp$aggregate(fact, fun, narm, opt)
		return (messages(x, "aggregate"))
	} else {
		out <- rast(x)
		nl <- nlyr(out)
		opt <- spatOptions()
		out@cpp <- out@cpp$aggregate(fact, "sum", TRUE, opt)
		out <- messages(out, "aggregate")
		dims <- x@cpp$get_aggregate_dims(fact)

		vtest <- values(x, dataframe=TRUE, row=1, nrows=dims[1], col=1, ncols=dims[2])
		vtest <- as.list(vtest)
		test <- sapply(vtest, fun)
		dm <- dim(test)
		do_transpose = FALSE
		if (!is.null(dm)) {
			do_transpose = TRUE
		}
		if (inherits(test, "list")) {
			error("aggregate", "fun returns a list")
		}
		fun_ret <- 1
		if (length(test) > nl) {
			if ((length(test) %% nl) == 0) {
				fun_ret <- length(test) / nl
				nlyr(out) <- nlyr(x) * fun_ret
			} else {
				error("aggregate", "cannot use this function")
			}
		}

		b <- blocks(x, 4)

		nr <- max(1, floor(b$nrows[1] / fact[1])) * fact[1]
		nrs <- rep(nr, floor(nrow(x)/nr))
		d <- nrow(x) - sum(nrs)
		if (d > 0) nrs <- c(nrs, d)
		b$row <- c(0, cumsum(nrs))[1:length(nrs)] + 1
		b$nrows <- nrs
		b$n <- length(nrs)
		outnr <- ceiling(b$nrows / fact[1]);
		outrows  <- c(0, cumsum(outnr))[1:length(outnr)] + 1
		nc <- ncol(x)


		if (inherits(cores, "cluster")) {
			doPar <- TRUE
		} else if (cores > 1) {
			doPar <- TRUE
			cores <- parallel::makeCluster(cores)
			on.exit(parallel::stopCluster(cores))
			export_args(cores, ..., caller="aggregate")
		} else {
			doPar <- FALSE
		}

		mpl <- prod(dims[5:6]) * fun_ret
		readStart(x)
		on.exit(readStop(x))
		ignore <- writeStart(out, filename, overwrite, sources=sources(x), wopt=wopt)
		if (doPar) {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc)
				v <- x@cpp$get_aggregates(v, b$nrows[i], dims)
				v <- parallel::parSapply(cores, v, fun, ...)
				if (length(v) != outnr[i] * mpl) {
					error("aggregate", "this function does not return the correct number of values")
				}
				if (do_transpose) {
					v <- t(v)
				}
				writeValues(out, v, outrows[i], outnr[i])
			}
		} else {
			for (i in 1:b$n) {
				v <- readValues(x, b$row[i], b$nrows[i], 1, nc)
				v <- x@cpp$get_aggregates(v, b$nrows[i], dims)
				v <- sapply(v, fun, ...)
				if (length(v) != outnr[i] * mpl) {
					error("aggregate", "this function does not return the correct number of values")
				}
				if (do_transpose) {
					v <- t(v)
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


aggregate_attributes <- function(d, by, fun=NULL, count=TRUE, ...) {
	i <- sapply(d, is.numeric)
	i[colnames(d) %in% by] <- FALSE
	j <- 1:length(by)
	da <- db <- NULL
	if (!is.null(fun)) {
		if (any(i)) {
			if (is.character(fun)) {
				f <- match.fun(fun)
				da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], f, ...)
				names(da)[-j] <- paste0(fun, "_", names(da)[-j])
			} else {
				da <- aggregate(d[, i,drop=FALSE], d[, by, drop=FALSE], fun, ...)
				names(da)[-j] <- paste0("agg_", names(da)[-j])
			}
		} else {
			da <- unique(d[, by, drop=FALSE])
		}
	} else {
		i[] <- FALSE
	}
	i[colnames(d) %in% by] <- TRUE
	if (any(!i)) {
		db <- aggregate(d[, !i,drop=FALSE], d[, by, drop=FALSE], .agg_uf)
		#db <- db[, colSums(is.na(db)) < nrow(db), drop=FALSE]
		if (NCOL(da)>1) {
			da <- merge(da, db, by=by)
		} else {
			da <- db
		}
	}

	if (count) {
		dn <- aggregate(d[, by[1],drop=FALSE], d[, by, drop=FALSE], length)
		colnames(dn)[ncol(dn)] = "agg_n"
		if (NCOL(da) > 1) {
			if (nrow(dn) > 0) {
				dn <- merge(da, dn, by=by)
			} else {
				dn <- da
				dn$agg_n <- 1
			}
		}
		dn
	} else {
		da
	}
}


setMethod("aggregate", signature(x="SpatVector"),
	function(x, by=NULL, dissolve=TRUE, fun="mean", count=TRUE, ...) {
		if (inherits(by, "SpatVector")) {
			error("use 'zonal' to aggregate a SpatVector with a SpatVector")
		}
		if (is.null(by)) {
			x$aggregate_by_variable = 1;
			x@cpp <- x@cpp$aggregate("aggregate_by_variable", dissolve)
			x$aggregate_by_variable = NULL;
		} else {
			if (is.character(by)) {
				by <- unique(by)
				iby <- match(by, names(x))
				if (any(is.na(iby))) {
					bad <- paste(by[is.na(iby)], collapse=", ")
					error("aggregate", "invalid name(s) in by: ", bad)
				}
			} else if (is.numeric(by)) {
				by <- unique(by)
				iby <- round(by)
				if (any((iby < 1) | (iby > ncol(x)))) {
					bad <- iby[(iby < 1) | (iby > ncol(x))]
					error("aggregate", "invalid column number in by: ", bad)
				}
			} else {
				error("aggregate", "by should be character or numeric")
			}

			d <- values(x)
			mvars <- FALSE
			if (length(iby) > 1) {
				cvar <- apply(d[, iby], 1, function(i) paste(i, collapse="_"))
				by <- basename(tempfile())
				values(x) <- NULL
				x[[by]] <- cvar
				mvars <- TRUE
			} else {
				by <- names(x)[iby]
			}

			x@cpp <- x@cpp$aggregate(by, dissolve)
			messages(x)

			if (mvars) {
				d[[by]] <- cvar
				a <- aggregate_attributes(d, c(by, names(d)[iby]), fun=fun, count=count, ...)
			} else {
				a <- aggregate_attributes(d, names(d)[iby], fun=fun, count=count, ...)
			}

			if (any(is.na(d[[by]]))) {
				# because NaN and NA are dropped
				i <- nrow(a)+(1:2)
				a[i,] <- c(NA, NaN)
			}
			i <- match(x[[by,drop=TRUE]], a[[by]])
			i <- i[!is.na(i)]
			if (mvars) {
				a[[by]] <- NULL
			}
			values(x) <- a[i,,drop=FALSE]
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

