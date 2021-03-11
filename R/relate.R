
setMethod("relate", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, relation) {
		out <- x@ptr$relate_between(y@ptr, relation)
		x <- messages(x, "relate")
		out[out == 2] <- NA
		matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
	}
)


setMethod("relate", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y, relation, ...) {
		y <- as.polygons(y)
		relate(x, y, relation, ...)
	}
)

setMethod("relate", signature(x="SpatExtent", y="SpatVector"), 
	function(x, y, relation, ...) {
		x <- as.polygons(x)
		relate(x, y, relation, ...)
	}
)


setMethod("relate", signature(x="SpatExtent", y="SpatExtent"), 
	function(x, y, relation, ...) {
		x <- as.polygons(x)
		y <- as.polygons(y)
		relate(x, y, relation, ...)
	}
)


setMethod("relate", signature(x="SpatVector", y="missing"), 
	function(x, y, relation, pairs=FALSE, symmetrical=FALSE) {
		out <- x@ptr$relate_within(relation, symmetrical)
		x <- messages(x, "relate")
		out[out == 2] <- NA
		if (symmetrical) {
			class(out) <- "dist"
			attr(out, "Size") <- nrow(x)
			attr(out, "Diag") <- FALSE
			attr(out, "Upper") <- FALSE
		} else {
			out <- matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
		}	
		if (pairs) {
			out <- mat2wide(out, symmetrical)
		}
		out
	}
)


setMethod("adjacent", signature(x="SpatRaster"), 
	function(x, cells, directions="rook", include=FALSE) {
		v <- x@ptr$adjacent(cells-1, directions, include)
		messages(x, "adjacent")
		v <- do.call(rbind, v)
		rownames(v) <- cells
		return(v+1)
	}
)


setMethod("adjacent", signature(x="SpatVector"), 
	function(x, type="rook", pairs=TRUE, symmetrical=FALSE) {
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
		if (pairs) {
			mat2wide(a, symmetrical, 1)
		} else {
			a
		}
	}
)



setMethod("near", signature(x="SpatVector"), 
	function(x, distance=0, k=1, centroids=TRUE, symmetrical=TRUE) {
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
