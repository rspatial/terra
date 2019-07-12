
setMethod("rastApply", signature(x="SpatRaster"), 
function(x, indices, fun, ..., filename="", overwrite=FALSE, wopt=list()) {

	stopifnot(!any(is.na(indices)))
	if (!is.factor(indices)) {
		indices <- as.factor(indices)
	}
	nms <- as.character(indices)
	ind <- as.integer(indices)
	d <- unique(data.frame(nms, ind, stringsAsFactors=FALSE))
	uin <- d[,2]
	nms <- make.names(d[,1])

	txtfun <- .makeTextFun(match.fun(fun))
	if (class(txtfun) == "character") { 
		if (txtfun %in% c("max", "min", "mean", "prod", "sum", "any", "all")) {
			opt <- .runOptions(filename, overwrite, wopt)
			na.rm <- isTRUE(list(...)$na.rm)
			x@ptr <- x@ptr$apply(indices, txtfun, na.rm, nms, opt)	
			return(show_messages(x))
		}		
	}

	nl <- nlyr(x)
	ind <- rep_len(indices, nl)
	out <- rast(x)
	nlyr(out) <- length(uin)
	names(out) <- nms
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, ncol(out), TRUE)
		v <- lapply(uin, function(i, ...) apply(v[, ind==uin[i], drop=FALSE], 1, fun, ...))
		v <- do.call(cbind, v)
		writeValues(out, v, b$row[i])
	}
	writeStop(out)
	readStop(x)
	return(out)
}
)


