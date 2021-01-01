
setMethod("relate", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, relation) {
		out <- x@ptr$relate_between(y@ptr, relation)
		x <- messages(x, "relate")
		out[out == 2] <- NA
		matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
	}
)


setMethod("relate", signature(x="SpatVector", y="ANY"), 
	function(x, y, relation) {
		yy <- try(vect(y), silent=TRUE)
		if (!inherits(yy, "SpatVector")) {
			yy <- try(ext(y), silent=TRUE)
			if (!inherits(yy, "SpatExtent")) {
				stop("cannot use argument 'y'")
			}
			yy <- as.polygons(yy)
		}
		relate(x, yy, relation)
	}
)

setMethod("relate", signature(x="ANY", y="SpatVector"), 
	function(x, y, relation) {
		xx <- try(vect(x), silent=TRUE)
		if (!inherits(xx, "SpatVector")) {
			xx <- try(ext(x), silent=TRUE)
			if (!inherits(xx, "SpatExtent")) {
				stop("cannot use argument 'x'")
			}
			xx <- as.polygons(xx)
		}
		relate(xx, y, relation)
	}
)

setMethod("relate", signature(x="ANY", y="ANY"), 
	function(x, y, relation) {
		xx <- try(vect(x), silent=TRUE)
		if (!inherits(xx, "SpatVector")) {
			xx <- try(ext(x), silent=TRUE)
			if (!inherits(xx, "SpatExtent")) {
				stop("cannot use argument 'x'")
			}
			xx <- as.polygons(xx)
		}
		yy <- try(vect(y), silent=TRUE)
		if (!inherits(yy, "SpatVector")) {
			yy <- try(ext(y), silent=TRUE)
			if (!inherits(yy, "SpatExtent")) {
				stop("cannot use argument 'y'")
			}
			yy <- as.polygons(xx)
		}
		relate(xx, yy, relation)
	}
)



setMethod("relate", signature(x="SpatVector", y="missing"), 
	function(x, y, relation, pairs=FALSE, symmetrical=FALSE) {
		out <- x@ptr$relate_within(relation, symmetrical)
		x <- messages(x, "relate")
		out[out == 2] <- NA
		out <- matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
		if (pairs) {
			out <- mat2wide(out, symmetrical)
		}
		out
	}
)


setMethod("adjacent", signature(x="SpatRaster"), 
	function(x, cells, directions="rook", include=FALSE, ...) {
		v <- x@ptr$adjacent(cells-1, directions, include)
		messages(x, "adjacent")
		v <- do.call(rbind, v)
		rownames(v) <- cells
		return(v+1)
	}
)


setMethod("adjacent", signature(x="SpatVector"), 
	function(x, type="rook", symmetrical=FALSE, ...) {
		type <- match.arg(tolower(type), c("intersects", "touches", "queen", "rook"))
		stopifnot(geomtype(x) == "polygons")
		a <- x@ptr$relate_within(type, TRUE)
		x <- messages(x, "relate")
		a[a == 2] <- NA
		class(a) <- "dist"
		attr(a, "Size") <- nrow(x)
		attr(a, "Diag") <- FALSE
		attr(a, "Upper") <- FALSE
		a <- as.matrix(a)
		mat2wide(a, symmetrical, 1)
	}
)



setMethod("near", signature(x="SpatVector"), 
	function(x, distance=0, k=1, centroids=TRUE, symmetrical=TRUE, ...) {
		if ((geomtype(x) == "polygons") && centroids) {
			x <- centroids(x)
		}
		if (distance > 0) {
			d <- distance(x, pairs=TRUE, symmetrical=symmetrical)
			d[d[,3] <= distance, 1:2,drop=FALSE]		
		} else {
			k <- max(1, min(round(k), nrow(x)))
			d <- as.matrix(distance(x, pairs=FALSE))
			diag(d) <- NA
			t(apply(d, 1, function(i) order(i)[1:k]))
		}
	}
)
