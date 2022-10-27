
setMethod("tapp", signature(x="SpatRaster"),
function(x, index, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list()) {

	stopifnot(!any(is.na(index)))

	if ((length(index) == 1) && is.character(index)) {
		choices <- c("years", "months", "week", "days", "doy", "yearmonths")
		i <- pmatch(tolower(index), choices)
		if (is.na(i)) {
			error("tapp", paste("invalid time step. Use one of:", paste(choices, collapse=", ")))
		}
		if (!x@ptr$hasTime) {
			error("tapp", "x has no time data")
		}
		choice <- choices[i]
		if (choice == "doy") {
			index <- format(time(x, "days"), "%j")
		} else if (choice == "week") {
			index <- strftime(time(x, "days"), format = "%V")
		} else {
			index <- time(x, choice)
			if (choice == "yearmonths") {
				year <- trunc(index)
				month <- 12 * (index - year) + 1
				year <- formatC(year, width=4, flag = "0")
				month <- formatC(month, width=2, flag = "0")
				index <- paste0(year, month)
			} else {
				index <- as.character(index)
			}
		}
		#time(x) <- NULL
	}

	nl <- nlyr(x)
	if (length(index) > nl) {
		error("tapp", "length(index) > nlyr(x)")
	}
	index <- rep_len(index, nl)
	if (!is.factor(index)) {
		index <- factor(index, levels=unique(index))
	}
	nms <- as.character(index)
	ind <- as.integer(index)
	d <- unique(data.frame(nms, ind, stringsAsFactors=FALSE))
	uin <- d[,2]
	nms <- make.names(d[,1])
	nms <- nms[uin]
	txtfun <- .makeTextFun(fun)
	if (inherits(txtfun, "character")) {
		if (txtfun %in% .cpp_funs) {
			opt <- spatOptions(filename, overwrite, wopt=wopt)
			narm <- isTRUE(list(...)$na.rm)
			x@ptr <- x@ptr$apply(index, txtfun, narm, nms, opt)
			return(messages(x, "tapp"))
		}
	}
	fun <- match.fun(fun)

	readStart(x)
	on.exit(readStop(x), add=TRUE)

	testnc <- min(ncol(x), 11)
	v <- readValues(x, 1, 1, 1, testnc, TRUE)

	test <- apply(v, 1, FUN=fun, ...)
	transpose = FALSE
	nlout <- 1
	if (ncol(test) > 1) {
		if (ncol(test) == testnc) {
			transpose = TRUE
			nlout <- nrow(test)
			addnms <- rownames(test)
		} else {
			nlout <- ncol(test)
			addnms <- colnames(test)
		}
		nms <- paste(rep(nms, each=length(addnms)), rep(addnms, length(nms)), sep="_")
	}

	out <- rast(x)
	nlyr(out) <- nlout * length(uin)
	names(out) <- nms

	doclust <- FALSE
	if (inherits(cores, "cluster")) {
		doclust <- TRUE
	} else if (cores > 1) {
		doclust <- TRUE
		cores <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cores))
	}


	b <- writeStart(out, filename, overwrite, sources=sources(x), wopt=wopt)

	if (doclust) {
		export_args(cores, ..., caller="tapp")
		pfun <- function(x, ...) apply(x, 1, FUN=fun, ...)
		parallel::clusterExport(cores, "pfun", environment())
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, ncol(out), TRUE)
			v <- lapply(uin, function(i) v[, ind==i, drop=FALSE])
			v <- parallel::parLapply(cores, v, pfun, ...)
			if (transpose) {
				v <- t(do.call(rbind, v))
			} else {
				v <- do.call(cbind, v)
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	} else {
		for (i in 1:b$n) {
			v <- readValues(x, b$row[i], b$nrows[i], 1, ncol(out), TRUE)
			# like this, na.rm is not passed to FUN
			# v <- lapply(uin, function(j, ...) apply(v[, ind==uin[j], drop=FALSE], 1, FUN=fun, ...))
			# like this it works
			v <- lapply(uin, function(j) apply(v[, ind==uin[j], drop=FALSE], 1, FUN=fun, ...))
			if (transpose) {
				v <- t(do.call(rbind, v))
			} else {
				v <- do.call(cbind, v)
			}
			writeValues(out, v, b$row[i], b$nrows[i])
		}
	}
	out <- writeStop(out)
	return(out)
}
)


