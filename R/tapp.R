


setMethod("tapp", signature(x="SpatRaster"),
function(x, index, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list()) {

	stopifnot(!any(is.na(index)))

	prename <- ""
	out_time <- double()
	out_tstep <- ""
	out_tz <- "UTC"
	if (length(index) == 1) {
		if (is.character(index)) {
			choices <- c("years", "months", "week", "days", "doy", "yearmonths", "yearweeks", "7days", "10days", "15days")
			i <- pmatch(tolower(index), choices)
			if (is.na(i)) {
				error("tapp", paste("invalid time step. Use one of:", paste(choices, collapse=", ")))
			}
			if (!x@cpp$hasTime) {
				error("tapp", "x has no time data")
			}
			choice <- choices[i]
			if (choice == "doy") {
				# or POSIXlt$yday
				index <- format(time(x, "days"), "%j")
				prename <- "doy_"
			} else if (choice == "week") {
				index <- strftime(time(x, "days"), format = "%V")
				prename <- "week_"
			} else if (choice == "yearweeks") {
				index <- yearweek(time(x, "days"))
				prename <- "yw_"
			} else if (choice == "7days") {
				index <- as.integer(format(time(x, "days"), "%j"))
				index <- as.character((index-1) %/% 7 + 1)
				prename <- "d7_"
			} else if (choice == "10days") {
				index <- as.integer(format(time(x, "days"), "%j"))
				index <- as.character((index-1) %/% 10 + 1)
				prename <- "d10_"
			} else if (choice == "15days") {
				index <- as.integer(format(time(x, "days"), "%j"))
				index <- as.character((index-1) %/% 15 + 1)
				prename <- "d15_"
			} else {
				index <- time(x, choice)	
				out_time <- time_as_seconds(x)[!duplicated(index)]
				out_tstep <- choice
				out_tz <- attr(index, "tzone")
				if (is.null(out_tz)) out_tz = "UTC"
				if (choice == "yearmonths") {
					year <- trunc(index)
					month <- 12 * (index - year) + 1
					year <- formatC(year, width=4, flag = "0")
					month <- formatC(month, width=2, flag = "0")
					index <- paste0(year, month)
					prename <- "ym_"
				} else {
					index <- as.character(index)
					if (choice == "months") {
						prename <- "m_"
					} else if (choice == "days") {
						prename <- "d_"
					} else if (choice == "years") {
						prename <- "y_"
					} 
				}
			}
		} else if (is.function(index)) {
			index <- as.character(index(time(x)))
		} 
	}

	nl <- nlyr(x)
	if (length(index) > nl) {
		error("tapp", "length(index) > nlyr(x)")
	} else if (length(unique(index)) == 1) {
		warn("tapp", "it is not sensible to a single value as index (use app instead)")	
	}
	index <- rep_len(index, nl)
	if (!is.factor(index)) {
		index <- factor(index, levels=unique(index))
	}
	nms <- paste0(prename, as.character(index))
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
			x@cpp <- x@cpp$apply(index, txtfun, narm, nms, out_time, out_tstep, out_tz, opt)
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
	if (NCOL(test) > 1) {
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
	if (out_tstep != "") {
		time(out, out_tstep) <- out_time 
	}
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


