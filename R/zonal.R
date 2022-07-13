
setMethod("zonal", signature(x="SpatRaster", z="SpatRaster"),
	function(x, z, fun="mean", ..., as.raster=FALSE, filename="", wopt=list())  {
		if (nlyr(z) > 1) {
			z <- z[[1]]
		}
		zname <- names(z)
		txtfun <- .makeTextFun(fun)
		if (inherits(txtfun, "character") && (txtfun %in% c("max", "min", "mean", "sum", "notNA", "isNA"))) {
			na.rm <- isTRUE(list(...)$na.rm)
			old <- isTRUE(list(...)$old)
			opt <- spatOptions()
			if (old) { # for testing, to be removed
				sdf <- x@ptr$zonal_old(z@ptr, txtfun, na.rm, opt)
			} else {
				sdf <- x@ptr$zonal(z@ptr, txtfun, na.rm, opt)
			}
			messages(sdf, "zonal")
			out <- .getSpatDF(sdf)
		} else {
			fun <- match.fun(fun)
			nl <- nlyr(x)
			res <- list()
			vz <- values(z)
			nms <- names(x)
			for (i in 1:nl) {
				d <- stats::aggregate(values(x[[i]]), list(zone=vz), fun, ...)
				colnames(d)[2] <- nms[i]
				res[[i]] <- d
			}
			out <- res[[1]]
			if (nl > 1) {
				for (i in 2:nl) {
					out <- merge(out, res[[i]])
				}
			}
		}
		if (as.raster) {
			if (is.null(wopt$names)) {
				wopt$names <- names(x)
			}
			levels(z) <- NULL
			subst(z, out[,1], out[,-1], filename=filename, wopt=wopt)
		} else {
			if (is.factor(z)) {
				levs <- active_cats(z)[[1]]
				m <- match(out$zone, levs[,1])
				out$zone <- levs[m, 2]
			}
			colnames(out)[1] <- zname
			out
		}
	}
)


setMethod("zonal", signature(x="SpatVector", z="SpatVector"),
	function(x, z, fun=mean, ..., weighted=FALSE, as.polygons=FALSE)  {
		if (geomtype(z) != "polygons") {
			error("zonal", "x must be points, and z must be polygons")
		}
		if (nrow(x) == 0) {
			error("zonal", "x is empty")		
		}
		isn <- which(sapply(values(x[1,]), is.numeric))
		if (!any(isn)) {
			error("zonal", "x has no numeric variables (attributes) to aggregate")
		}
		x <- x[,isn]
		if (geomtype(x) == "points") {
			r <- !relate(x, z, "disjoint")
			i <- apply(r, 1, \(i) if(any(i)) which(i) else (NA))
			if (length(i) == 0) {
				error("zonal", "there are no points in x that overlap with the polygons in z")		
			}
			a <- aggregate(values(x), data.frame(zone=i), fun, ...)
		} else {
			if (as.polygons) {
				zz <- z
				values(zz) <- data.frame(zone = 1:nrow(zz))
				i <- intersect(zz, x)
			} else {
				values(z) <- data.frame(zone = 1:nrow(z))
				i <- intersect(z, x)			
			}
			if (nrow(i) == 0) {
				error("zonal", "the intersection of x and z is empty")			
			}
			v <- values(i)
			if (weighted) {
				if (geomtype(i) == "lines") {
					v$w <- perim(i)
				} else {
					v$w <- expanse(i)
				}
				s <- split(v, v$zone)
				n <- ncol(v)-2
				s <- lapply(s, function(d) {
						out <- rep(NA, n)
						for (i in 2:n) {
							out[i-1] <- weighted.mean(d[[i]], w = d$w)
						}
						out
					})
				a <- data.frame(as.integer(names(s)), do.call(rbind, s))
				colnames(a) <- names(v)[-ncol(v)]
			} else {
				a <- aggregate(v[,-1,drop=FALSE], v[,1,drop=FALSE], fun, ...)
			}
		}
		if (as.polygons) {
			f <- basename(tempfile())
			z[[f]] <- 1:nrow(z)
			names(a)[1] = f
			a <- merge(z, a, by=f, all.x=TRUE)
			a[[f]] <- NULL
		}
		a
	}
)


setMethod("global", signature(x="SpatRaster"),
	function(x, fun="mean", weights=NULL, ...)  {

		nms <- names(x)
		nms <- make.unique(nms)
		txtfun <- .makeTextFun(fun)

		opt <- spatOptions()
		if (!is.null(weights)) {
			stopifnot(inherits(weights, "SpatRaster"))
			stopifnot(txtfun %in% c("mean", "sum"))
			na.rm <- isTRUE(list(...)$na.rm)
			ptr <- x@ptr$global_weighted_mean(weights@ptr, txtfun, na.rm, opt)
			messages(ptr, "global")
			res <- (.getSpatDF(ptr))
			rownames(res) <- nms
			return(res)
		}

		if (inherits(txtfun, "character")) {
			if (txtfun %in% c("prod", "max", "min", "mean", "sum", "range", "rms", "sd", "sdpop", "notNA", "isNA")) {
				na.rm <- isTRUE(list(...)$na.rm)
				ptr <- x@ptr$global(txtfun, na.rm, opt)
				messages(ptr, "global")
				res <- .getSpatDF(ptr)

				rownames(res) <- nms
				return(res)
			}
		}

		nl <- nlyr(x)
		res <- list()
		for (i in 1:nl) {
			res[[i]] <- fun(values(x[[i]]), ...)
		}
		res <- do.call(rbind, res)
		res <- data.frame(res)

		# more efficient but more risky:
		#apply(data.frame(x), 2, fun, ...)

		if ((ncol(res) == 1) && (colnames(res) == "res")) {
			colnames(res) <- "global"
		}

		rownames(res) <- nms
		res
	}
)
