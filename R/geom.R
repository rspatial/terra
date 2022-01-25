

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



setMethod("split", signature(x="SpatVector"), 
	function(x, f) {
		if (length(f) > 1) {
			x <- copy(x)
			x$f <- f
			f <- "f"
		}
		x <- messages(x@ptr$split(f), "split")
		as.list.svc(x)
	}
)


setMethod("cover", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, identity=FALSE) {
		x@ptr <- x@ptr$cover(y@ptr, identity[1])
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
		x@ptr <- x@ptr$erase(y@ptr)
		messages(x, "erase")
	}
)

setMethod("erase", signature(x="SpatVector", y="missing"), 
	function(x) {
		x@ptr <- x@ptr$erase_self()
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
		x@ptr <- x@ptr$intersect(y@ptr)
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
		x@ptr <- x@ptr$crop_ext(y@ptr)
		x
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
		messages(x, "intersect")
	}
)


#setMethod("intersect", signature(x="SpatRaster", y="SpatRaster"),
#	function(x, y) {
#		a <- crop(x, y)
#		b <- crop(y, x)
#		c(a, b)
#	}
#)

setMethod("buffer", signature(x="SpatVector"), 
	function(x, width, quadsegs=10) {
		x@ptr <- x@ptr$buffer(width, quadsegs)
		messages(x, "buffer")
	}
)


setMethod("crop", signature(x="SpatVector", y="ANY"), 
	function(x, y) {
		if (inherits(y, "SpatVector")) {
			if (length(y) > 1) {
				y <- aggregate(y)
			}
			x@ptr <- x@ptr$crop_vct(y@ptr)
		} else {
			if (!inherits(y, "SpatExtent")) {
				y <- try(ext(y), silent=TRUE)
				if (inherits(y, "try-error")) {
					stop("y does not have a SpatExtent")
				}
			}
			x@ptr <- x@ptr$crop_ext(y@ptr)
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
	function(x) {
		x@ptr <- x@ptr$disaggregate()
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
		x@ptr <- x@ptr$rotate(angle, x0[1], y0[1])
		messages(x, "spin")
	}
)


setMethod("delauny", signature(x="SpatVector"), 
	function(x, tolerance=0, as.lines=FALSE) {
		x@ptr <- x@ptr$delauny(tolerance, as.lines)
		messages(x, "delauny")
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
	function(x, tolerance=0.1) {
		preserveTopology <- TRUE
		x@ptr <- x@ptr$simplify(tolerance, preserveTopology)
		messages(x, "simplifyGeom")
	}
)


setMethod("sharedPaths", signature(x="SpatVector"), 
	function(x) {
		x@ptr <- x@ptr$shared_paths()
		messages(x, "sharedPaths")
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
