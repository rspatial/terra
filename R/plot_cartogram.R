
setMethod("cartogram", signature(x="SpatVector"),
	function(x, var, type="nc")  {
		if (geomtype(x) != "polygons") {
			error("cartogram", "x must be polygons")
		}
		type <- match.arg(tolower(type), "nc")
		stopifnot(var %in% names(x))
		v <- as.numeric(as.vector(x[[var, drop=TRUE]]))
		i <- !is.na(v)
		x <- x[i]
		v <- v[i]

		i <- v > 0
		x <- x[i]
		v <- v[i]
		
		if (nrow(x) == 0) return(vect("POLYGON EMPTY"))

		f <- v / max(v)
		cxy <- crds(centroids(x, inside=TRUE))
		r <- lapply(1:length(v), function(i) {
			rescale(x[i,], f[i], x0=cxy[i,1], y0=cxy[i,2])
		})
		do.call(rbind, r)
	}
)

