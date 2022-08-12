
setMethod("rapp", signature(x="SpatRaster"),
function(x, first, last, fun, ..., allyrs=FALSE, fill=NA, clamp=FALSE, circular=FALSE, filename="", overwrite=FALSE, wopt=list()) {

	stopifnot(hasValues(x))
	firstval <- lastval <- NA
	if (inherits(first, "SpatRaster")) {
		first <- first[[1]]
		stopifnot(hasValues(first))
	} else {
		if (!is.numeric(first)) {
			error("rapp", "argument `first` should be numeric or SpatRaster")
		}
		firstval <- first
		stopifnot(first %in% 1:nlyr(x))
	}
	if (inherits(last, "SpatRaster")) {
		last <- last[[1]]
		stopifnot(hasValues(last))
	} else {
		if (!is.numeric(last)) {
			error("rapp", "argument `last` should be numeric or SpatRaster")
		}
		lastval <- last
		stopifnot(last %in% 1:nlyr(x))
	}
	if (!(is.na(firstval)) && (!(is.na(lastval)))) {
		error("rapp", "argument `first` or `last` must be a SpatRaster. Or use `app`")
	}
	if (!is.na(firstval)) {
		index <- last;
	} else if (!is.na(lastval)) {
		index <- first
	} else {
		index <- c(first, last)
	}
	compareGeom(x, index, lyrs=FALSE, crs=FALSE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE)

	if (!allyrs) {
		txtfun <- .makeTextFun(match.fun(fun))
		if (inherits(txtfun, "character")) {
			if (txtfun %in% .cpp_funs) {
				opt <- spatOptions(filename, overwrite, wopt=wopt)
				na.rm <- isTRUE(list(...)$na.rm)
				x@ptr <- x@ptr$rapply(index@ptr, firstval, lastval, txtfun, clamp, na.rm, circular, opt)
				return(messages(x, "rapp"))
			}
		}
	}
	out <- rast(x)
	v <- x@ptr$rappvals(index@ptr, firstval, lastval, clamp, allyrs, fill, 0, 1, circular)
	v <- sapply(v, fun, ...)
	if (is.list(v)) { error("rapp", "values returned by 'fun' do not have the same length for each cell") }
	nc <- ncol(out)
	trans = FALSE
	if (NCOL(v) == nc) {
		trans = TRUE
		nlyr(out) <- 	nrow(v)
	} else if (NROW(v) == nc) {
		nlyr(out) <- NCOL(v)
	} else if (length(v) == nc) {
		nlyr(out) <- 1
	}
	b <- writeStart(out, filename, overwrite, sources=sources(x), wopt=wopt, n=nlyr(x)*3)
	for (i in 1:b$n) {
		v <- x@ptr$rappvals(index@ptr, firstval, lastval, clamp, allyrs, fill, b$row[i]-1, b$nrows[i], circular)
		v <- sapply(v, fun, ...)
		if (trans) v = t(v)
		writeValues(out, as.vector(v), b$row[i], b$nrows[i])
	}
	out <- writeStop(out)
	return(out)
}
)


