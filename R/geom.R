

roundtrip <- function(x, coll=FALSE) {
	if (coll) {
		p <- methods::new("SpatVectorCollection")
		p@ptr <- x@ptr$bienvenue()
		return(p)
	} else {
		x@ptr <- x@ptr$allerretour()
		return(x)
	}
}

get_invalid_coords <- function(x) {
	x <- x[!x[,1], ]
	if (nrow(x) > 0) {
		id <- as.integer(rownames(x))
		txt <- x[,2]
		txt <- gsub("Ring Self-intersection\\[", "", txt)
		txt <- gsub("Self-intersection\\[", "", txt)
		txt <- gsub("Too few points in geometry component\\[", "", txt)
		txt <- unlist(strsplit(gsub("]", "", txt), " "))
		txt <- matrix(as.numeric(txt), ncol=2, byrow=TRUE)
		v <- vect(txt)
		values(v) <- data.frame(id=id, msg=x[,2])
		v
	} else {
		vect()
	}
}

setMethod("is.valid", signature(x="SpatVector"),
	function(x, messages=FALSE, as.points=FALSE) {
		if (as.points) messages = TRUE
		if (messages) {
			r <- x@ptr$geos_isvalid_msg()
			d <- data.frame(matrix(r, ncol=2, byrow=TRUE))
			d[,1] = d[,1] == "\001"
			colnames(d) <- c("valid", "reason")
			if (as.points) {
				p <- try(get_invalid_coords(d), silent=TRUE)
				if (inherits(p, "try-error")) {
					warn("is.valid", "as.points failed, returning matrix")
					return(d)
				} else {
					return(p)
				}
			}
			d
		} else {
			x@ptr$geos_isvalid()
		}
	}
)

setMethod("makeValid", signature(x="SpatVector"),
	function(x) {
		x@ptr <- x@ptr$make_valid2()
		messages(x)
	}
)


setMethod("na.omit", signature("SpatVector"),
	function(object, field=NA, geom=FALSE) {
		if (geom) {
			g <- geom(object)
			g <- g[is.na(g[,"x"]) | is.na(g[,"y"]), 1]
			if (length(g) > 0) {
				object <- object[-g, ]
			}
		}
		if (!is.na(field)) {
			v <- values(object)
			if (field != "") {
				v <- v[, field, drop=FALSE]
			}
			i <- apply(v, 1, function(i) any(is.na(i)))
			if (any(i)) {
				object <- object[!i, ]
			}
		}
		object
	}
)


setMethod("deepcopy", signature("SpatVector"),
	function(x) {
		x@ptr <- x@ptr$deepcopy()
		x
	}
)


as.list.svc <- function(x) {
	v <- vect()
	lapply(1:x$size(),
		function(i) {
			v@ptr <- x$get(i-1)
			v
		})
}



setMethod("cover", signature(x="SpatVector", y="SpatVector"),
	function(x, y, identity=FALSE, expand=TRUE) {
		x@ptr <- x@ptr$cover(y@ptr, identity[1], expand[1])
		messages(x, "cover")
	}
)


setMethod("symdif", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		x@ptr <- x@ptr$symdif(y@ptr)
		messages(x, "symdif")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		x@ptr <- x@ptr$erase_agg(y@ptr)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="missing"),
	function(x, sequential=TRUE) {
		x@ptr <- x@ptr$erase_self(sequential)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatExtent"),
	function(x, y) {
		y <- as.polygons(y)
		x@ptr <- x@ptr$erase(y@ptr)
		messages(x, "erase")
	}
)

setMethod("gaps", signature(x="SpatVector"),
	function(x) {
		x@ptr <- x@ptr$gaps()
		messages(x, "gaps")
	}
)


setMethod("union", signature(x="SpatVector", y="missing"),
	function(x, y) {
		x@ptr <- x@ptr$union_self()
		messages(x, "union")
	}
)


setMethod("union", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		x@ptr <- x@ptr$union(y@ptr)
		messages(x, "union")
	}
)

setMethod("union", signature(x="SpatVector", y="SpatExtent"),
	function(x, y) {
		y <- as.vector(y)
		x@ptr <- x@ptr$union(y@ptr)
		messages(x, "union")
	}
)

setMethod("union", signature(x="SpatExtent", y="SpatExtent"),
	function(x, y) {
		x + y
	}
)


setMethod("intersect", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		x@ptr <- x@ptr$intersect(y@ptr, TRUE)
		messages(x, "intersect")
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatExtent"),
	function(x, y) {
		x@ptr = x@ptr$intersect(y@ptr)
		if (!x@ptr$valid) {
			return(NULL)
		}
		x
	}
)

setMethod("intersect", signature(x="SpatVector", y="SpatExtent"),
	function(x, y) {
		#x@ptr <- x@ptr$crop_ext(y@ptr)
		#x
		crop(x, y)
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatVector"),
	function(x, y) {
		y <- ext(y)
		x * y
	}
)

setMethod("mask", signature(x="SpatVector", mask="SpatVector"),
	function(x, mask, inverse=FALSE) {
		x@ptr <- x@ptr$mask(mask@ptr, inverse)
		messages(x, "mask")
	}
)


#setMethod("mask", signature(x="SpatVector", mask="sf"),
#	function(x, mask, inverse=FALSE) {
#		mask(x, vect(mask), inverse=inverse)
#	}
#)

#setMethod("intersect", signature(x="SpatRaster", y="SpatRaster"),
#	function(x, y) {
#		a <- crop(x, y)
#		b <- crop(y, x)
#		c(a, b)
#	}
#)

setMethod("buffer", signature(x="SpatVector"),
	function(x, width, quadsegs=10, capstyle="round", joinstyle="round", mitrelimit=NA, singlesided=FALSE) {
		if (is.character(width)) {
			if (!(width %in% names(x))) {
				error("buffer", paste(width, "is not a field in x"))
			}
			width <- x[[width, drop=TRUE]]
		}
		if (!is.numeric(width)) {
			error("buffer", "width is not numeric")
		}
		x@ptr <- x@ptr$buffer(width, quadsegs, tolower(capstyle), tolower(joinstyle), mitrelimit, singlesided)
		messages(x, "buffer")
	}
)


setMethod("crop", signature(x="SpatVector", y="ANY"),
	function(x, y, ext=FALSE) {
		if (ext) {
			y <- ext(y)
			x@ptr <- x@ptr$crop_ext(y@ptr, TRUE)
			return(x)
		}
		if (inherits(y, "SpatVector")) {
			x@ptr <- x@ptr$crop_vct(y@ptr)
		} else {
			if (!inherits(y, "SpatExtent")) {
				y <- try(ext(y), silent=TRUE)
				if (inherits(y, "try-error")) {
					stop("y does not have a SpatExtent")
				}
			}
			## crop_ext does not include points on the borders
			## https://github.com/rspatial/raster/issues/283
			#x@ptr <- x@ptr$crop_ext(y@ptr)
			y <- as.polygons(y)
			x@ptr <- x@ptr$crop_vct(y@ptr)
		}
		messages(x, "crop")
	}
)


setMethod("convHull", signature(x="SpatVector"),
	function(x, by="") {
		x@ptr <- x@ptr$hull("convex", by[1])
		messages(x, "convHull")
	}
)

setMethod("minRect", signature(x="SpatVector"),
	function(x, by="") {
		x@ptr <- x@ptr$hull("minrot", by[1])
		messages(x, "minRect")
	}
)


setMethod("disagg", signature(x="SpatVector"),
	function(x, segments=FALSE) {
		x@ptr <- x@ptr$disaggregate(segments[1])
		messages(x, "disagg")
	}
)




setMethod("flip", signature(x="SpatVector"),
	function(x, direction="vertical") {
		d <- match.arg(direction, c("vertical", "horizontal"))
		x@ptr <- x@ptr$flip(d == "vertical")
		messages(x, "flip")
	}
)



setMethod("spin", signature(x="SpatVector"),
	function(x, angle, x0, y0) {
		e <- as.vector(ext(x))
		if (missing(x0)) {
			x0 <- mean(e[1:2])
		}
		if (missing(y0)) {
			y0 <- mean(e[3:4])
		}
		angle <- angle[1]
		stopifnot(is.numeric(angle) && !is.nan(angle))
		x@ptr <- x@ptr$rotate(angle, x0, y0)
		messages(x, "spin")
	}
)


setMethod("delaunay", signature(x="SpatVector"),
	function(x, tolerance=0, as.lines=FALSE) {
		x@ptr <- x@ptr$delaunay(tolerance, as.lines)
		messages(x, "delaunay")
	}
)



voronoi_deldir <- function(x, bnd=NULL, eps=1e-09, ...){

	xy <- crds(x)
	dat <- values(x)
	if (nrow(dat > 0)) {
		dups <- duplicated(xy)
		if (any(dups)) {
			xy <- xy[!dups, ,drop=FALSE]
			dat <- dat[!dups, ,drop=FALSE]
		}
	} else {
		xy <- stats::na.omit(xy[, 1:2])
		xy <- unique(xy)
	}

	e <- bnd
	if (!is.null(e)) {
		e <- as.vector(ext(bnd))
	}

	dd <- deldir::deldir(xy[,1], xy[,2], rw=e, eps=eps, suppressMsge=TRUE)
	g <- lapply(deldir::tile.list(dd), function(i) cbind(i$ptNum, 1, i$x, i$y))
	g <- do.call(rbind, g)
	g <- vect(g, "polygons", crs=crs(x))
	if (nrow(g) == nrow(dat)) {
		values(g) <- dat
	} else {
		values(g) <- data.frame(id=dd$ind.orig)
	}
	g
}



setMethod("voronoi", signature(x="SpatVector"),
	function(x, bnd=NULL, tolerance=0, as.lines=FALSE, deldir=FALSE) {
		if (geomtype(x) != "points") {
			x <- as.points(x)
		}
		if (deldir) {
			voronoi_deldir(x, bnd, tolerance=tolerance)
		} else {
			if (is.null(bnd)) {
				bnd <- vect()
			} else {
				bnd <- as.polygons(ext(bnd))
			}
			x@ptr <- x@ptr$voronoi(bnd@ptr, tolerance, as.lines)
			messages(x, "voronoi")
		}
	}
)

setMethod("elongate", signature(x="SpatVector"),
	function(x, length=1, flat=FALSE) {
		x@ptr <- x@ptr$elongate(length, flat)
		messages(x, "elongate")
	}
)


setMethod("width", signature(x="SpatVector"),
	function(x, as.lines=FALSE) {
		x@ptr <- x@ptr$width()
		x <- messages(x, "width")
		if (!as.lines) {
			x <- perim(x)
		}
		x
	}
)


setMethod("clearance", signature(x="SpatVector"),
	function(x, as.lines=FALSE) {
		x@ptr <- x@ptr$clearance()
		x <- messages(x, "clearance")
		if (!as.lines) {
			x <- perim(x)
		}
		x
	}
)




setMethod("mergeLines", signature(x="SpatVector"),
	function(x) {
		x@ptr <- x@ptr$line_merge()
		messages(x, "line_merge")
	}
)

setMethod("makeNodes", signature(x="SpatVector"),
	function(x) {
		x@ptr <- x@ptr$make_nodes()
		messages(x, "makeNodes")
	}
)

setMethod("removeDupNodes", signature(x="SpatVector"),
	function(x, digits=-1) {
		x@ptr <- x@ptr$remove_duplicate_nodes(digits)
		messages(x, "removeDupNodes")
	}
)


setMethod("simplifyGeom", signature(x="SpatVector"),
	function(x, tolerance=0.1, preserveTopology=TRUE, makeValid=TRUE) {
		x@ptr <- x@ptr$simplify(tolerance, preserveTopology)
		x <- messages(x, "simplifyGeom")
		if (makeValid) {
			x <- makeValid(x)
		}
		x
	}
)

setMethod("thinGeom", signature(x="SpatVector"),
	function(x, threshold=1e-6, makeValid=TRUE) {
		x@ptr <- x@ptr$thin(threshold)
		x <- messages(x, "thinGeom")
		if (makeValid) {
			x <- makeValid(x)
		}
		x
	}
)

setMethod("sharedPaths", signature(x="SpatVector"),
	function(x, y=NULL) {
		if (is.null(y)) {
			x@ptr <- x@ptr$shared_paths(TRUE)
		} else {
			x@ptr <- x@ptr$shared_paths2(y@ptr, TRUE)
		}
		x <- messages(x, "sharedPaths")
		# sort data to ensure consistent order with spatial indices
		if (nrow(x) > 0) x <- x[order(x$id1, x$id2), ]
		x
	}
)


setMethod("snap", signature(x="SpatVector"),
	function(x, y=NULL, tolerance) {
		if (is.null(y)) {
			x@ptr <- x@ptr$snap(tolerance)
		} else {
			x@ptr <- x@ptr$snapto(y@ptr, tolerance)
		}
		messages(x, "snap")
	}
)


setMethod("combineGeoms", signature(x="SpatVector", y="SpatVector"),
	function(x, y, overlap=TRUE, boundary=TRUE, distance=TRUE, append=TRUE, minover=0.1, maxdist=Inf, dissolve=TRUE, erase=TRUE) {

		if ((geomtype(x) != "polygons") || (geomtype(y) != "polygons")) {
			error("combineGeoms", "x and y must be polygons")
		}
		if (nrow(x) == 0) {
			if (append) {
				return(rbind(x, y))
			} else {
				return(x)
			}
		}
		if (nrow(y) == 0) {
			return(x)
		}

		xcrs <- crs(x)
		ycrs <- crs(y)
		if ((xcrs == "") || (ycrs == "")) {
			error("combineGeoms", "x and y must have a crs")
		} else if (xcrs != ycrs) {
			error("combineGeoms", "x and y do not have the same crs")
		}

		dx <- values(x)
		dy <- values(y)
		values(x) = data.frame(idx=1:nrow(x))
		values(y) = data.frame(idy=1:nrow(y))
		y <- erase(y) # no self-overlaps
		if (overlap) {
			#avoid Warning message: [intersect] no intersection
 			xy <- suppressWarnings(intersect(y, x))
			if (nrow(xy) > 0) {
				xy$aint <- expanse(xy)
				a <- values(xy)
				a <- a[order(a$idy, -a$aint),]
				a <- a[!duplicated(a$idy),]
				yi <- y[a$idy,]
				atot <- expanse(yi)
				a <- a[(a$aint / atot) >= minover, ]
				if (nrow(a) > 0) {
					if (erase) {
						ye <- erase(y, x)
						i <- na.omit(match(a$idy, ye$idy))
						if (length(i) > 0) {
							yi <- ye[i,]
							values(yi) <- data.frame(idx=a$idx[i])
						} else {
							yi <- vect()
						}
					} else {
						yi <- y[a$idy,]
						values(yi) <- data.frame(idx=a$idx)
					}
					if (nrow(yi) > 0) {
						x <- aggregate(rbind(x, yi), "idx", dissolve=dissolve, counts=FALSE)
					}
					y <- y[-a$idy,]
				}
			}
		}

		if (boundary && (nrow(y) > 0)) {
			ye <- erase(y, x)
			p <- sharedPaths(ye, x)
			if (nrow(p) > 0) {
				p$s <- perim(p)
				p <- values(p)
				p <- p[order(p$id1, -p$s),]
				p <- p[!duplicated(p$id1),]
				if (erase) {
					i <- p$id1
					yi <- ye[p$id1,]
				} else {
					i <- ye$idy[p$id1]
					i <- match(i, y$idy)
					yi <- y[i,]
				}
				yi$idx <- 0
				yi$idx[i] <- p$id2
				yi$idy <- NULL
				x <- aggregate(rbind(x, yi), "idx", dissolve=dissolve, counts=FALSE)
				y <- y[-i,]
			}
		}

		if (distance && (nrow(y) > 0) && (maxdist > 0)) {
			n <- nearest(y, x)
			n <- n[n$distance <= maxdist, ]
			if (nrow(n) > 0) {
				yi <- y[n$from_id, ]
				yi$idx <- n$to_id
				yi$idy <- NULL
				x <- aggregate(rbind(x, yi), "idx", dissolve=FALSE, counts=FALSE)
				y <- y[-n$from_id, ]
			}
		}

		values(x) <- dx[x$idx, ,drop=FALSE]
		if (append && (nrow(y) > 0)) {
			values(y) <- dy[y$idy, ,drop=FALSE]
			if (erase) {
				y <- erase(y, x)
			}
			x <- rbind(x, y)
		}
		x
	}
)



setMethod("split", signature(x="SpatVector", f="ANY"),
	function(x, f) {
		if (length(f) > 1) {
			x <- deepcopy(x)
			x$f <- f
			f <- "f"
		}
		x <- messages(x@ptr$split(f), "split")
		as.list.svc(x)
	}
)


setMethod("split", signature(x="SpatVector", f="SpatVector"),
	function(x, f) {
		if (geomtype(x) != "polygons") error("split", "first argument must be polygons")
		if (geomtype(f) != "lines") error("split", "second argument must be lines")
		values(f) <- NULL
		u <- union(ext(f), ext(x)) + 10
		e <- ext(u)
		md <- max(e$xmax-e$xmin, e$ymax-e$ymin)
		lin <- elongate(f, md)
		xf <- rbind(as.lines(x), f)
		xf <- aggregate(xf)
		m <- makeNodes(xf)
		p <- as.polygons(m)
		intersect(x, p)
	}
)
