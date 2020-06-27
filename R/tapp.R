
setMethod("tapp", signature(x="SpatRaster"), 
function(x, index, fun, ..., filename="", overwrite=FALSE, wopt=list()) {

	stopifnot(!any(is.na(index)))
	if (!is.factor(index)) {
		index <- as.factor(index)
	}
	nms <- as.character(index)
	ind <- as.integer(index)
	d <- unique(data.frame(nms, ind, stringsAsFactors=FALSE))
	uin <- d[,2]
	nms <- make.names(d[,1])

	txtfun <- .makeTextFun(fun)
	if (inherits(txtfun, "character")) { 
		if (txtfun %in% c("max", "min", "mean", "prod", "sum", "any", "all")) {
			opt <- .runOptions(filename, overwrite, wopt)
			narm <- isTRUE(list(...)$na.rm)
			x@ptr <- x@ptr$apply(index, txtfun, narm, nms, opt)	
			return(show_messages(x, "tapp"))
		}		
	}
	fun <- match.fun(fun)

	nl <- nlyr(x)
	ind <- rep_len(index, nl)
	out <- rast(x)
	nlyr(out) <- length(uin)
	names(out) <- nms
	stopifnot(readStart(x))
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, ncol(out), TRUE)
		# like this, na.rm is not passed to FUN
		# v <- lapply(uin, function(j, ...) apply(v[, ind==uin[j], drop=FALSE], 1, FUN=fun, ...))
		# like this it works
		v <- lapply(uin, function(j) apply(v[, ind==uin[j], drop=FALSE], 1, FUN=fun, ...))
		v <- do.call(cbind, v)
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	readStop(x)
	out <- writeStop(out)
	return(out)
}
)


