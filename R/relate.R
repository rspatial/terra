
setMethod("is.related", signature(x="SpatVector", y="SpatVector"),
	function(x, y, relation) {
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)

setMethod("is.related", signature(x="SpatVector", y="SpatExtent"),
	function(x, y, relation) {
		y <- as.polygons(y)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)

setMethod("is.related", signature(x="SpatExtent", y="SpatVector"),
	function(x, y, relation) {
		x <- as.polygons(x)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)


setMethod("is.related", signature(x="SpatVector", y="SpatRaster"),
	function(x, y, relation) {
		y <- as.polygons(y, ext=TRUE)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)

setMethod("is.related", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, relation) {
		x <- as.polygons(x, ext=TRUE)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)


setMethod("is.related", signature(x="SpatExtent", y="SpatRaster"),
	function(x, y, relation) {
		x <- as.polygons(x)
		y <- as.polygons(y, ext=TRUE)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)

setMethod("is.related", signature(x="SpatRaster", y="SpatExtent"),
	function(x, y, relation) {
		x <- as.polygons(x, ext=TRUE)
		y <- as.polygons(y)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)

setMethod("is.related", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, relation) {
		x <- as.polygons(x, ext=TRUE)
		y <- as.polygons(y, ext=TRUE)
		out <- x@pntr$is_related(y@pntr, relation)
		x <- messages(x, "is.related")
		out
	}
)


setMethod("relate", signature(x="SpatVector", y="SpatVector"),
	function(x, y, relation, pairs=FALSE, na.rm=TRUE) {
		if (pairs) {
			out <- x@pntr$related_between(y@pntr, relation[1], na.rm[1])
			messages(x, "relate")
			if (length(out[[1]]) == 0) {
				cbind(id.x=0,id.y=0)[0,,drop=FALSE]
			} else {
				names(out) <- c("id.x", "id.y")
				do.call(cbind, out) + 1
			}
		} else {
			out <- x@pntr$related_between(y@pntr, relation[1], TRUE)
			messages(x, "relate")
			m <- matrix(FALSE, nrow(x), nrow(y))
			if (length(out[[1]]) > 0) {
				m[do.call(cbind, out) + 1] <- TRUE
			}
			m

#			out <- x@pntr$relate_between(y@pntr, relation, TRUE, TRUE)
#			messages(x, "relate")
#			out[out == 2] <- NA
#			matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
		}
	}
)

#setMethod("which.related", signature(x="SpatVector", y="SpatVector"),
#	function(x, y, relation) {
#		out <- x@pntr$which_related(y@pntr, relation)
#		x <- messages(x, "which.related")
#		out <- do.call(cbind, out) + 1
#		colnames(out) <- c("id.x", "id.y")
#		out
#	}
#)



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


setMethod("relate", signature(x="SpatVector", y="SpatRaster"),
	function(x, y, relation, ...) {
		y <- as.polygons(y, ext=TRUE)
		relate(x, y, relation, ...)
	}
)

setMethod("relate", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, relation, ...) {
		x <- as.polygons(x, ext=TRUE)
		relate(x, y, relation, ...)
	}
)

setMethod("relate", signature(x="SpatExtent", y="SpatRaster"),
	function(x, y, relation, ...) {
		x <- as.polygons(x)
		y <- as.polygons(y, ext=TRUE)
		relate(x, y, relation, ...)
	}
)

setMethod("relate", signature(x="SpatRaster", y="SpatExtent"),
	function(x, y, relation, ...) {
		x <- as.polygons(x, ext=TRUE)
		y <- as.polygons(y)
		relate(x, y, relation, ...)
	}
)

setMethod("relate", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, relation, ...) {
		x <- as.polygons(x, ext=TRUE)
		y <- as.polygons(y, ext=TRUE)
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
	function(x, y, relation, pairs=FALSE, na.rm=TRUE) {

		if (pairs) {
			out <- x@pntr$related_within(relation, na.rm[1])
			messages(x, "relate")
			if (length(out[[1]]) == 0) {
				cbind(id.1=0,id.2=0)[0,,drop=FALSE]
			} else {
				names(out) <- c("id.1", "id.2")
				do.call(cbind, out) + 1
			}
		} else {
			out <- x@pntr$related_within(relation, TRUE)
			messages(x, "relate")
			out <- do.call(cbind, out) + 1
			m <- matrix(FALSE, nrow(x), nrow(x))
			m[out] <- TRUE
			m
			#if (symmetrical) {
			#	as.dist(m)
			#} else {
			#	m
			#}
		}

		#out <- x@pntr$relate_within(relation, symmetrical)
		#out[out == 2] <- NA
		#if (symmetrical) {
		#	class(out) <- "dist"
		#	attr(out, "Size") <- nrow(x)
		#	attr(out, "Diag") <- FALSE
		#	attr(out, "Upper") <- FALSE
		#} else {
		#	out <- matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
		#}
		#if (pairs) {
		#	out <- mat2wide(out, symmetrical)
		#}
		#out
	}
)



setMethod("adjacent", signature(x="SpatRaster"),
	function(x, cells, directions="rook", pairs=FALSE, include=FALSE, symmetrical=FALSE) {
		cells <- cells - 1
		if (inherits(directions, "matrix")) {
			directions[!is.finite(directions)] <- 0
			if (isTRUE(all(directions == 0))) {
				error("adjacent", "directions are all FALSE")
			}
			v <- x@pntr$adjacentMat(cells, as.logical(t(directions)), dim(directions), include)
		} else {
			#if (pairs) include <- FALSE
			v <- x@pntr$adjacent(cells,  as.character(directions)[1], include)
		}
		messages(x, "adjacent")
		if (pairs) {
			v <- cbind(from=rep(cells, each=length(v)/length(cells)), to=v)
			v <- v[!is.na(v[,2]), ]
			if (symmetrical) {
				#v <- unique(cbind(pmin(v[,1], v[,2]), pmax(v[,1], v[,2])))
				v <- .unique_symmetric_rows(v[,1], v[,2])
			}
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
		a <- x@pntr$relate_within(type, TRUE)
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
	function(x, y=NULL, distance=0, k=1, centroids=TRUE, symmetrical=TRUE, method="geo") {
		
		k <- round(k)
		if (distance <= 0 && k < 1) {
			error("nearby", "either distance or k must be a positive number")
		}

		if (!(method %in% c("geo", "haversine", "cosine"))) {
			error("nearby", "not a valid method. Should be one of: 'geo', 'haversine', 'cosine'")
		}

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
				d <- distance(x, y, method=method)
				d <- cbind(from_id=rep(1:nrow(d), ncol(d)), to_id=rep(1:ncol(d), each=nrow(d)), distance=as.vector(d))
			} else {
				d <- distance(x, pairs=TRUE, symmetrical=symmetrical, method=method)
			}
			d[d[,3] <= distance, 1:2, drop=FALSE]
		} else {
			if (hasy) {
				k <- max(1, min(round(k), (nrow(y)-1)))
			} else {
				k <- max(1, min(round(k), (nrow(x)-1)))
			}
#			if (k > 1) {
				if (hasy) {
					d <- distance(x, y)
				} else {
					d <- as.matrix(distance(x, pairs=FALSE, method=method))
					diag(d) <- NA
				}
				d <- t(apply(d, 1, function(i) order(i)[1:k]))
				if (k==1) d <- t(d)
				d <- cbind(1:length(x), d)
#			} else {
#				d <- nearest(x)
#				d <- values(d)[, c("from_id", "to_id")]
#			}
			colnames(d) <- c("id", paste0("k", 1:k))
			d
		}
	}
)



setMethod("nearest", signature(x="SpatVector"),
	function(x, y=NULL, pairs=FALSE, centroids=TRUE, lines=FALSE, method="geo") {
		if ((geomtype(x) == "polygons") && centroids) {
			x <- centroids(x)
		}
		within <- FALSE
		if (is.null(y)) {
			within <- TRUE
			y <- x
		} else if ((geomtype(y) == "polygons") && centroids) {
			y <- centroids(y)
		}
		z <- x
		if (!(method %in% c("geo", "haversine", "cosine"))) {
			error("nearest", "not a valid method. Should be one of: 'geo', 'haversine', 'cosine'")
		}
		
		if (within) {
			z@pntr <- x@pntr$near_within(method)
		} else {
			z@pntr <- x@pntr$near_between(y@pntr, pairs, method)
		}
		z <- messages(z, "nearest")
		if (geomtype(z) == "points") { 
		
			if (within) {
				names(z)[1] <- "to_id"
				z$to_id <- z$to_id + 1
				crd <- crds(x)
				X <- cbind(crd[,1], crd[z$to_id, 1])
				Y <- cbind(crd[,2], crd[z$to_id, 2])
				if (lines) {
					geom <- cbind(rep(1:nrow(x), each=2), 1, as.vector(t(X)), as.vector(t(Y)))
					z <- vect(geom, "lines", crs=crs(z), atts=values(z))
				}
				values(z) <- data.frame(from_id=1:nrow(z), from_x=X[,1], from_y=Y[,2], to_id=z$to_id, to_x=X[,2], to_y=Y[,2], distance=z$distance)
				return(z)
			}
		
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
				z$to_id = z$to_id + 1
				z$from_id = z$from_id + 1
				return(z)
			}
		} else {
			if (lines) return(z)
			values(y) <- data.frame(to_id=1:nrow(y))
			dis <- perim(z)
			zz <- as.points(z)
			from <- zz[seq(1, nrow(zz), 2), ]
			to <- zz[seq(2, nrow(zz), 2), ]
			values(to) <- data.frame(id=1:nrow(to))
			to_int <- as.data.frame(intersect(to, y))
			to_int <- to_int[!duplicated(to_int[,1]), ]
			to_int <- as.data.frame(merge(to, to_int, by="id", all.x=TRUE))
			if (any(is.na(to_int$to_id))) {
				zz <- as.points(elongate(z, 1))
				from2 <- zz[seq(1, nrow(zz), 2), ]
				to2 <- zz[seq(2, nrow(zz), 2), ]
				values(to2) <- data.frame(id=1:nrow(to2))
				to_int2 <- as.data.frame(intersect(to2, y))
				colnames(to_int2)[2] <- "to_id2"
				to_int <- merge(to_int, to_int2, all.x=TRUE)
				i <- is.na(to_int$to_id)
				to_int$to_id[i] <- to_int$to_id2[i]
				to_int <- to_int[,1:2]
			}
			
			to_int <- to_int[order(to_int[["id"]]), ]

			if (nrow(to_int) > nrow(to)) {
				to_int <- aggregate(to_int[, "to_id",drop=FALSE], to_int[,"id",drop=FALSE], function(x)x[1])
#			}
#			if (nrow(to_int) < nrow(to)) {
#				to_int <- rep(NA, nrow(to))
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

