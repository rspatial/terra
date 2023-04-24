
setMethod("cartogram", signature(x="SpatVector"),
	function(x, var, type)  {
		if (geomtype(x) != "polygons") {
			error("cartogram", "x must be polygons")
		}
		type <- match.arg(tolower(type), "nc")
		stopifnot(var %in% names(x))
		v <- as.numeric(as.vector(x[[var, drop=TRUE]]))
		if (!any(!is.na(v))) stop(paste("no numeric values in", var))
		if (any(v <= 0)) stop(paste("non-positive values in", var))
		x <- x[!is.na(v)]
		v <- v[!is.na(v)]
		f <- v / max(v)
		cxy <- crds(centroids(x, inside=TRUE))
		r <- lapply(1:length(v), function(i) {
			rescale(x[i,], f[i], x0=cxy[i,1], y0=cxy[i,2])
		})
		do.call(rbind, r)
	}
)

