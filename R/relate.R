
setMethod("is.related", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, relation) {
		out <- x@ptr$is_related(y@ptr, relation)
		x <- messages(x, "is.related")
		out
	}
)

setMethod("is.related", signature(x="SpatVector", y="SpatExtent"), 
	function(x, y, relation) {
		y <- as.polygons(y)
		out <- x@ptr$is_related(y@ptr, relation)
		x <- messages(x, "is.related")
		out
	}
)



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
	function(x, cells, directions="rook", pairs=FALSE, include=FALSE) {
		cells <- cells - 1
		if (inherits(directions, "matrix")) {
			v <- x@ptr$adjacentMat(cells, as.logical(directions), dim(directions), include)
		} else {
			if (pairs) include <- FALSE
			v <- x@ptr$adjacent(cells,  as.character(directions)[1], include)
		}
		messages(x, "adjacent")
		if (pairs) {
			v <- cbind(from=rep(cells, each=length(v)/length(cells)), to=v)
			v <- v[!is.na(v[,2]), ]
		} else {
			v <- matrix(v, nrow=length(cells), byrow=TRUE)
			if (!include) rownames(v) <- cells
		}
		v + 1
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


setMethod("nearby", signature(x="SpatVector"), 
	function(x, y=NULL, distance=0, k=1, centroids=TRUE, symmetrical=TRUE) {
		if ((geomtype(x) == "polygons") && centroids) {
			x <- centroids(x)
		}
		hasy <- !is.null(y)
		if (hasy) {
			if ((geomtype(y) == "polygons") && centroids) {
				y <- centroids(y)
			}
		}
		if (distance > 0) {
			if (hasy) {
				d <- distance(x, y)
				d <- cbind(from_id=rep(1:nrow(d), ncol(d)), to_id=rep(1:ncol(d), each=nrow(d)), distance=as.vector(d))
			} else {
				d <- distance(x, pairs=TRUE, symmetrical=symmetrical)
			}
			d[d[,3] <= distance, 1:2, drop=FALSE]
		} else {
			if (hasy) {
				k <- max(1, min(round(k), (nrow(y)-1)))
			} else {
				k <- max(1, min(round(k), (nrow(x)-1)))
			}
			if (k > 1) {
				if (hasy) {
					d <- distance(x, y)
				} else {
					d <- as.matrix(distance(x, pairs=FALSE))
					diag(d) <- NA
				}
				d <- t(apply(d, 1, function(i) order(i)[1:k]))
				if (k==1) d <- t(d)
				d <- cbind(1:length(x), d)
			} else {
				d <- nearest(x)
				d <- values(d)[, c("from_id", "to_id")]
			}
			colnames(d) <- c("id", paste0("k", 1:k))
			d
		}
	}
)



setMethod("nearest", signature(x="SpatVector"), 
	function(x, y=NULL, pairs=FALSE, centroids=TRUE, lines=FALSE) {
		if ((geomtype(x) == "polygons") && centroids) {
			x <- centroids(x)
		}
		within <- FALSE
		if (is.null(y)) {
			within <- TRUE
			y <- x
		} else {
			if ((geomtype(y) == "polygons") && centroids) {
				y <- centroids(y)
			}
		}
		z <- x
		if (within) {
			z@ptr <- x@ptr$near_within()
		} else {
			z@ptr <- x@ptr$near_between(y@ptr, pairs)
		}
		z <- messages(z, "nearest")
		if (geomtype(z) == "points") { #lonlat points
			if (lines) {
				x <- z[,c(2,5), drop=TRUE]
				y <- z[,c(3,6), drop=TRUE]
				geom <- cbind(rep(1:nrow(x), each=2), 1, as.vector(t(x)), as.vector(t(y)))
				zz <- vect(geom, "lines", crs=crs(z))
				values(zz) <- values(z)
				zz$to_id = zz$to_id + 1
				zz$from_id = zz$from_id + 1
				return(zz)
			} else {
				values(z) <- data.frame(from_id = 1:nrow(z), to_id = z$id + 1, distance=z$distance)
				return(z)
			}
		} else {
			if (lines) return(z)
			dis <- perim(z)
			z <- as.points(z)
			from <- z[seq(1, nrow(z), 2), ]
			to <- z[seq(2, nrow(z), 2), ]
			values(to) <- data.frame(id=1:nrow(to))
			values(y) <- data.frame(to_id=1:nrow(y))
			to_int <- as.data.frame(intersect(to, y))
			to_int <- to_int[order(to_int[["id"]]), ]
			if (nrow(to_int) > nrow(to)) {
				to_int <- aggregate(to_int[, "to_id",drop=FALSE], to_int[,"id",drop=FALSE], function(x)x[1]) 
			} 
			if (nrow(to_int) < nrow(to)) {
				to_int <- rep(NA, nrow(to))
			} else {
				to_int <- to_int[,2] 
			}
			from <- geom(from)[, c("x", "y"),drop=FALSE]
			to <- geom(to)[, c("x", "y"),drop=FALSE]
			d <- data.frame(1:nrow(from), from, to_int, to, dis)
			colnames(d) <- c("from_id", "from_x", "from_y", "to_id", "to_x", "to_y", "distance")
			vect(d, c("to_x", "to_y"), crs=crs(x))
		}
	}
)

