
setMethod("cartogram", signature(x="SpatVector"),
	function(x, var, type="nc", inside=FALSE, exp=1)  {
		if (geomtype(x) != "polygons") {
			error("cartogram", "x must be polygons")
		}
		type <- tolower(type[1])
		type <- match.arg(type, c("nc", "circles"))
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
		if (exp <= 0) exp <- 1
		exp <- sqrt(exp)
		cntrds <- centroids(x, inside=inside)
		if (type == "circles") {
			w <- max(width(buffer(cntrds, width(x))))
			ff <- exp * (sqrt(f / pi) * w/4)
			buffer(cntrds, ff)
		} else {
			cxy <- crds(cntrds)
			s <- sqrt(expanse(x))
			ff <- s[which.max(f)] / s
			ff <- ff * sqrt(f) * exp
			r <- lapply(1:length(v), function(i) {
				rescale(x[i,], ff[i], x0=cxy[i,1], y0=cxy[i,2])
			})
			b = do.call(rbind, r)
		}
	}
)

