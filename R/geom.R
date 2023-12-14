
buffer2 <- function(x, width, quadsegs=10) {
	if (is.character(width)) {
		if (!(width %in% names(x))) {
			error("buffer2", paste(width, "is not a field in x"))
		}
		width <- x[[width, drop=TRUE]]
	}
	if (!is.numeric(width)) {
		error("buffer2", "width is not numeric")
	}
	x@cpp <- x@cpp$buffer2(width, quadsegs)
	messages(x, "buffer2")
}


buffer3 <- function(x, width, quadsegs=10) {
	if (is.character(width)) {
		if (!(width %in% names(x))) {
			error("buffer3", paste(width, "is not a field in x"))
		}
		width <- x[[width, drop=TRUE]]
	}
	if (!is.numeric(width)) {
		error("buffer3", "width is not numeric")
	}
	x@cpp <- x@cpp$buffer3(width, quadsegs)
	messages(x, "buffer3")
}

buffer4 <- function(x, width, quadsegs=10) {
	if (is.character(width)) {
		if (!(width %in% names(x))) {
			error("buffer4", paste(width, "is not a field in x"))
		}
		width <- x[[width, drop=TRUE]]
	}
	if (!is.numeric(width)) {
		error("buffer4", "width is not numeric")
	}
	x@cpp <- x@cpp$buffer4(width, quadsegs)
	messages(x, "buffer4")
}

roundtrip <- function(x, coll=FALSE) {
	if (coll) {
		p <- methods::new("SpatVectorCollection")
		p@cpp <- x@cpp$bienvenue()
		return(p)
	} else {
		x@cpp <- x@cpp$allerretour()
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
			r <- x@cpp$geos_isvalid_msg()
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
			x@cpp$geos_isvalid()
		}
	}
)

setMethod("makeValid", signature(x="SpatVector"),
	function(x) {
		x@cpp <- x@cpp$make_valid2()
		messages(x)
	}
)



setMethod("is.na", signature(x="SpatVector"),
	function(x) {
		x@cpp$naGeoms()
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
		stopifnot(is.vector(field))
		if (!is.na(field[1])) {
			field <- field[!is.na(field)]
			v <- values(object)
			if (!any(field == "")) {
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
		x@cpp <- x@cpp$deepcopy()
		x
	}
)


as.list.svc <- function(x) {
	v <- vect()
	lapply(1:x$size(),
		function(i) {
			v@cpp <- x$get(i-1)
			v
		})
}



setMethod("cover", signature(x="SpatVector", y="SpatVector"),
	function(x, y, identity=FALSE, expand=TRUE) {
		x@cpp <- x@cpp$cover(y@cpp, identity[1], expand[1])
		messages(x, "cover")
	}
)


setMethod("symdif", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		x@cpp <- x@cpp$symdif(y@cpp)
		messages(x, "symdif")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		x@cpp <- x@cpp$erase_agg(y@cpp)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="missing"),
	function(x, sequential=TRUE) {
		x@cpp <- x@cpp$erase_self(sequential)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="SpatExtent"),
	function(x, y) {
		y <- as.polygons(y)
		x@cpp <- x@cpp$erase(y@cpp)
		messages(x, "erase")
	}
)

setMethod("gaps", signature(x="SpatVector"),
	function(x) {
		x@cpp <- x@cpp$gaps()
		messages(x, "gaps")
	}
)


setMethod("union", signature(x="SpatVector", y="missing"),
	function(x, y) {
		x@cpp <- x@cpp$union_self()
		messages(x, "union")
	}
)


setMethod("union", signature(x="SpatVector", y="SpatVector"),
	function(x, y) {
		if (geomtype(x) != "polygons") {
			unique(rbind(x, y))
		} else {
			x@cpp <- x@cpp$union(y@cpp)
			messages(x, "union")
		}
	}
)

setMethod("union", signature(x="SpatVector", y="SpatExtent"),
	function(x, y) {
		y <- as.vector(y)
		x@cpp <- x@cpp$union(y@cpp)
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
		x@cpp <- x@cpp$intersect(y@cpp, TRUE)
		messages(x, "intersect")
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatExtent"),
	function(x, y) {
		x@cpp <- x@cpp$intersect(y@cpp)
		if (!x@cpp$valid_notempty) {
			return(NULL)
		}
		x
	}
)

setMethod("intersect", signature(x="SpatVector", y="SpatExtent"),
	function(x, y) {
		#x@cpp <- x@cpp$crop_ext(y@cpp)
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
		x@cpp <- x@cpp$mask(mask@cpp, inverse)
		messages(x, "mask")
	}
)

setMethod("intersect", signature(x="SpatExtent", y="SpatRaster"),
	function(x, y) {
		x <- align(x, y, snap="near")
		intersect(x, ext(y))
	}
)

setMethod("intersect", signature(x="SpatRaster", y="SpatExtent"),
	function(x, y) {
		intersect(y, x)
	}
)




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
		x@cpp <- x@cpp$buffer(width, quadsegs, tolower(capstyle), tolower(joinstyle), mitrelimit, singlesided)
		messages(x, "buffer")
	}
)



setMethod("crop", signature(x="SpatVector", y="ANY"),
	function(x, y, ext=FALSE) {
		if (ext) {
			y <- ext(y)
			x@cpp <- x@cpp$crop_ext(y@cpp, TRUE)
			return(x)
		}
		if (inherits(y, "SpatVector")) {
			x@cpp <- x@cpp$crop_vct(y@cpp)
		} else {
			if (!inherits(y, "SpatExtent")) {
				y <- try(ext(y), silent=TRUE)
				if (inherits(y, "try-error")) {
					stop("y does not have a SpatExtent")
				}
			}
			## crop_ext does not include points on the borders
			## https://github.com/rspatial/raster/issues/283
			#x@cpp <- x@cpp$crop_ext(y@cpp)
			y <- as.polygons(y)
			x@cpp <- x@cpp$crop_vct(y@cpp)
		}
		messages(x, "crop")
	}
)


setMethod("convHull", signature(x="SpatVector"),
	function(x, by="") {
		x@cpp <- x@cpp$hull("convex", by[1])
		messages(x, "convHull")
	}
)

setMethod("minRect", signature(x="SpatVector"),
	function(x, by="") {
		x@cpp <- x@cpp$hull("minrot", by[1])
		messages(x, "minRect")
	}
)

setMethod("minCircle", signature(x="SpatVector"),
	function(x, by="") {
		x@cpp <- x@cpp$hull("circle", by[1])
		messages(x, "minCircle")
	}
)


setMethod("disagg", signature(x="SpatVector"),
	function(x, segments=FALSE) {
		x@cpp <- x@cpp$disaggregate(segments[1])
		messages(x, "disagg")
	}
)




setMethod("flip", signature(x="SpatVector"),
	function(x, direction="vertical") {
		d <- match.arg(direction, c("vertical", "horizontal"))
		x@cpp <- x@cpp$flip(d == "vertical")
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
		x@cpp <- x@cpp$rotate(angle, x0, y0)
		messages(x, "spin")
	}
)


setMethod("delaunay", signature(x="SpatVector"),
	function(x, tolerance=0, as.lines=FALSE) {
		x@cpp <- x@cpp$delaunay(tolerance, as.lines)
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
		if (nrow(x) ==0) {
			error("voronoi", "input has no geometries")
		}
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
			x@cpp <- x@cpp$voronoi(bnd@cpp, tolerance, as.lines)
			messages(x, "voronoi")
		}
	}
)

setMethod("elongate", signature(x="SpatVector"),
	function(x, length=1, flat=FALSE) {
		x@cpp <- x@cpp$elongate(length, flat)
		messages(x, "elongate")
	}
)


setMethod("width", signature(x="SpatVector"),
	function(x, as.lines=FALSE) {
		x@cpp <- x@cpp$width()
		x <- messages(x, "width")
		if (!as.lines) {
			x <- perim(x)
		}
		x
	}
)


setMethod("clearance", signature(x="SpatVector"),
	function(x, as.lines=FALSE) {
		x@cpp <- x@cpp$clearance()
		x <- messages(x, "clearance")
		if (!as.lines) {
			x <- perim(x)
		}
		x
	}
)




setMethod("mergeLines", signature(x="SpatVector"),
	function(x) {
		x@cpp <- x@cpp$line_merge()
		messages(x, "line_merge")
	}
)

setMethod("makeNodes", signature(x="SpatVector"),
	function(x) {
		x@cpp <- x@cpp$make_nodes()
		messages(x, "makeNodes")
	}
)

setMethod("removeDupNodes", signature(x="SpatVector"),
	function(x, digits=-1) {
		x@cpp <- x@cpp$remove_duplicate_nodes(digits)
		messages(x, "removeDupNodes")
	}
)


setMethod("simplifyGeom", signature(x="SpatVector"),
	function(x, tolerance=0.1, preserveTopology=TRUE, makeValid=TRUE) {
		x@cpp <- x@cpp$simplify(tolerance, preserveTopology)
		x <- messages(x, "simplifyGeom")
		if (makeValid) {
			x <- makeValid(x)
		}
		x
	}
)

setMethod("thinGeom", signature(x="SpatVector"),
	function(x, threshold=1e-6, makeValid=TRUE) {
		x@cpp <- x@cpp$thin(threshold)
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
			x@cpp <- x@cpp$shared_paths(TRUE)
		} else {
			x@cpp <- x@cpp$shared_paths2(y@cpp, TRUE)
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
			x@cpp <- x@cpp$snap(tolerance)
		} else {
			x@cpp <- x@cpp$snapto(y@cpp, tolerance)
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
						i <- stats::na.omit(match(a$idy, ye$idy))
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
					yi$idx <- p$id2
					yi$idy <- NULL
				} else {
					i <- ye$idy[p$id1]
					i <- match(i, y$idy)
					yi <- y[i,]
					yi$idx <- 0
					yi$idx[i] <- p$id2[i]
				}
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
		x <- messages(x@cpp$split(f), "split")
		as.list.svc(x)
	}
)


setMethod("split", signature(x="SpatVector", f="SpatVector"),
	function(x, f) {
		if (geomtype(x) != "polygons") error("split", "first argument must be polygons")
		if (!(geomtype(f) %in% c("lines", "polygons"))) {
			error("split", "argument 'f' must have a lines or polygons geometry")
		}
		values(f) <- NULL
		ex <- ext(x)
		i <- intersect(ext(x), ext(f))
		if (is.null(i)) {
			warn("split", "x and f do not intersect")
			return(x)
		}
		r <- relate(x, f, "intersects")
		if (sum(r) == 0) {
			warn("split", "x and f do not intersect")
			return(x)
		}
		r <- r[rowSums(r) > 0, ,drop=FALSE]
		y <- x	
		values(y) <- NULL
		ss <- vector("list", nrow(r))
		if (geomtype(f) == "lines") {
			for (i in 1:nrow(r)) {
				yi <- y[i]
				yi <- disagg(yi)
				add <- NULL
				if (nrow(yi) > 1) {
					ri <- relate(yi, f, "intersects")
					if (any(!ri)) {
						add <- yi[!ri]
						yi <- yi[ri]
					}
				}
				lin <- intersect(f[r[i,],], yi)
				v <- rbind(as.lines(yi), lin)
				v <- makeNodes( aggregate(v) )
				v <- as.polygons(v)
				if (!is.null(add)) {
					v <- combineGeoms(v, add, overlap=FALSE, boundary=FALSE, distance=TRUE, dissolve=FALSE, erase=FALSE)
				}
				v$id <- i
				ss[[i]] <- v
			}
			ss <- vect(ss)
			v <- values(x)
			values(ss) <- v[ss$id, ]
			rbind(x[-i,], ss)
		} else { #if (geomtype(f) == "polygons") {
			for (i in 1:nrow(r)) {
				yi <- y[i]
				v <- rbind(intersect(yi, f[r[i,],]),
							   erase(yi, f[r[i,],]))
				v$id <- i
				ss[[i]] <- v
			}
			ss <- vect(ss)
			v <- values(x)
			values(ss) <- v[ss$id, ]
			rbind(x[-i,], ss)
		
		}
	}
)



setMethod("forceCCW", signature(x="SpatVector"),
	function(x) {
		x <- deepcopy(x)
		x@cpp$make_CCW()
		messages(x)
	}
)


setMethod("is.empty", signature(x="SpatVector"),
	function(x) {
		nrow(x) == 0
	}
)
