
setMethod("cartogram", signature(x="SpatVector"), 
	function(x, var, type)  {
		type <- match.arg(tolower(type), "nc")
		stopifnot(var %in% names(x))
		v <- as.numeric(as.vector(x[[var, drop=TRUE]]))
		if (!any(!is.na(v))) stop(paste("no numeric values in", var))
		if (any(v <= 0)) stop(paste("non-positive values in", var))
		x <- x[!is.na(v)]
		v <- v[!is.na(v)]
		f <- v / max(v)
		r <- lapply(1:length(v), function(i) rescale(x[i,], f[i]))
		do.call(c, r)
	}
)

