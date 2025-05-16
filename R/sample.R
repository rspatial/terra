

.km_regular <- function(x, n, ...) {
	x <- x[[1]]
	crs(x) <- "local" # for the search_radius
	z <- init(x, "xy")
	if (hasValues(x)) {
		z <- mask(z, x)
	}
	km <- k_means(z, n, ...)
	kpols <- as.polygons(km)
	a <- centroids(kpols, inside=FALSE)
	cds <- crds(a)
	cells <- cellFromXY(x, cds)
	i <- is.na(x[cells])
	if (any(i)) {
		j <- which(i)
		opt <- spatOptions()
		zz <- deepcopy(x)
		crs(zz) <- crs(kpols) <- "local" 
		for (k in j) {
			y <- crop(zz, kpols[k], mask=TRUE)
			ye <- as.vector(ext(y))
			sr <- max(diff(ye[1:2]), diff(ye[3:4]))
			e <- y@pntr$extractBuffer(cds[k,1], cds[k,2], sr, opt)
			if (!is.na(e[[1]])) {
				cells[k] <- cellFromXY(x, xyFromCell(y, e[[3]] + 1))
			} else {
				b <- centroids(kpols[k], inside=TRUE)
				bcds <- crds(b)
				cells[k] <- cellFromXY(x, bcds)
			}
		}
	}
	cells
}


regular_exact <- function(r, size) {
	size <- round(size)
	stopifnot(size > 0)
	if (size >= ncell(r)) {
		return(xyFromCell(r, 1:ncell(r)))
	}
	x <- (ymax(r) - ymin(r)) / (xmax(r) - xmin(r))
	if (x < 1) {
		nc <- min(max(1, size * 1/(1 + 1/x)), ncol(r))
		nr <- (size / nc)
	} else {
		nr <- min(max(1, size * 1/(1 + x)), nrow(r))
		nc <- (size / nr)
	}

	addrow <- addcol <- 0
	floor <- TRUE

	if ((nc%%1) >= (nr%%1)) {
		nc <- round(nc)
		if ((nc * ceiling(nr)) <= size) {
			nr <- ceiling(nr)
			floor <- FALSE
		} else {
			nr <- floor(nr)
			if ((nr * nc) > size) {
				nc <- nc-1
			}
		}
		d <- size - nr * nc
		if (d > 0) {
			if (floor) {
				addcol <- d
			} else {
				addrow <- d 
			}
		}
	} else {
		nr <- round(nr)
		if ((nr * ceiling(nc)) <= size) {
			nc <- ceiling(nc)
			floor <- FALSE
		} else {
			nc <- floor(nc)
			if ((nr * nc) > size) {
				nr <- nr-1
			}
		}
		d <- size - nr * nc
		if (d > 0) {
			if (floor) {
				addcol <- d
			} else {
				addrow <- d	
			}
		}
	}

	if (addcol > 0) {
		ncs <- rep(nr, nc+1)
		b <- (nr+addcol)/2
		ncs[1] <- ceiling(b)
		ncs[nc+1] <- floor(b)
		nc=ncs
	} else if (addrow > 0) {
		nrs <- rep(nc, nr+1)
		b <- (nc+addrow)/2
		nrs[1] <- ceiling(b)
		nrs[nr+1] <- floor(b)
		nr=nrs
	}
	nnr <- length(nr)
	nnc <- length(nc)
	
	if ((nnr == 1) && (nnc == 1)) {
		dx <- (xmax(r) - xmin(r)) / nc
		x <- xmin(r) + dx/2 + (0:(nc-1)) * dx 
		dy <- (ymax(r) - ymin(r)) / nr
		y <- ymin(r) + dy/2 + (0:(nr-1)) * dy 
		expand.grid(x=x, y=y)
	} else if (nnc == 1) {	
		dy <- (ymax(r) - ymin(r)) / nnr
		y <- ymin(r) + dy/2 + (0:(nnr-1)) * dy 
		out <- vector("list", nnr)
		for (i in 1:nnr) {
			dx <- (xmax(r) - xmin(r)) / nr[i]
			x <- xmin(r) + dx/2 + (0:(nr[i]-1)) * dx 
			out[[i]] <- expand.grid(x=x, y=y[i])
		}
		do.call(rbind, out)
	} else {
		dx <- (xmax(r) - xmin(r)) / nnc
		x <- xmin(r) + dx/2 + (0:(nnc-1)) * dx 
		out <- vector("list", nnc)
		for (i in 1:nnc) {
			dy <- (ymax(r) - ymin(r)) / nc[i]
			y <- ymin(r) + dy/2 + (0:(nc[i]-1)) * dy 
			out[[i]] <- expand.grid(x=x[i], y=y)
		}
		do.call(rbind, out)
	}
}



sampleWeights <- function(x, size, replace=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL) {

	if (nlyr(x) > 1) {
		x <- x[[1]]
		warn("spatSample", "only the first layer of x is used")
	}
	if (!hasValues(x)) {
		error("spatSample", "x has no values")
	}

	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	x <- classify(x, cbind(-Inf, 0, NA))
	res <- as.data.frame(x, cells=cells, xy=(xy | as.points))
	if (!replace) {
		if (size >= nrow(res)) {
			i <- 1:nrow(res)
		} else {
			i <- sample.int(nrow(res), size, prob=res[,ncol(res)], replace=replace)
		}
	} else {
		i <- sample.int(nrow(res), size, prob=res[,ncol(res)], replace=replace)
	}
	res <- res[i, , drop=FALSE]
	if (!values) res <- res[ , 1:(cells + 2*(xy | as.points)), drop=FALSE]
	if (as.points) {
		res <- vect(res, c("x", "y"), crs=crs(x))
		if (!xy) {
			res$x <- NULL
			res$y <- NULL
		}
	} #else if (as.df) {
		#res <- data.frame(res)
	#}
	res
}


sampleStratMemory <- function(x, size, replace, lonlat, ext=NULL, weights=NULL, warn=TRUE, each) {
	if (!is.null(ext)) {
		xold <- rast(x)
		x <- crop(x, ext)
		cells <- cells(xold, ext)
		if (!is.null(weights)) {
			weights <- crop(weights, ext)
		}
	} else {
		cells <- 1:ncell(x)	
	}
	
	doprob <- TRUE
	prob <- NULL
	if (!is.null(weights)) {
		if (!inherits(weights, "SpatRaster")) {
			error("spatSample", "weights must be a SpatRaster")
		}
		if (!compareGeom(x, weights)) {
			error("spatSample", "geometry of weights does not match the geometry of x")
		}
		v <- na.omit(cbind(cell=cells, values(x), values(weights)))		
	} else if (lonlat) {
		v <- cbind(cell=cells, values(x), abs(cos(pi * values(init(x, "y")) / 360)))
	} else {
		v <- cbind(cell=cells, values(x))
		doprob <- FALSE
	}

	v <- v[!is.na(v[,2]), ]
	uv <- sort(unique(v[,2]))
	nuv <- length(uv)
	if (each) {
		sz <- rep(size, nuv)
	} else {	
		sz <- rep(floor(size / nuv), nuv)
		d <- size - sum(sz)
		i <- sample(nuv, d)
		sz[i] <- sz[i] + 1
	}
	ys <- vector(mode="list", length=nuv)
		
	for (i in seq_len(length(uv))) {
		if (sz[i] == 0) next
		vv <- v[v[,2] == uv[i], ,drop=FALSE]
		if (doprob) prob <- vv[,3]
		if (replace) {
			s <- sample.int(nrow(vv), sz[i], prob=prob, replace=TRUE)
		} else {
			s <- sample.int(nrow(vv), min(sz[i], nrow(vv)), prob=prob, replace=FALSE)
		}
		ys[[i]] <- vv[s,-3]
	}
	ys <- do.call(rbind, ys)

	if (warn) {
		ta <- tapply(ys[,1], ys[,2], length)
		sz <- sz[sz > 0]
		ta <- names(ta)[ta < sz]
		if (length(ta) > 0) {
			warn("spatSample", 'fewer samples than requested are available for group(s): ', paste(ta, collapse=', '))
		}
	}

	ys
}



sampleStratified_old <- function(x, size, replace=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=TRUE, xy=FALSE, ext=NULL, warn=TRUE, exp=5, weights=NULL, exhaustive=FALSE, lonlat, each) {

	if (nlyr(x) > 1) {
		x <- x[[1]]
		warn("spatSample", "only the first layer of x is used")
	}
	if (!hasValues(x)) {
		error("spatSample", "x has no values")
	}


	if ((blocks(x, n=4)$n == 1) || exhaustive) {
		res <- sampleStratMemory(x, size, replace, lonlat, ext, weights, warn, each=each)
	} else {
		f <- unique(x)[,1]
		exp <- max(1, exp)
		ss <- exp * size * length(f)

		if (is.null(weights)) {
			if ((!lonlat) && (ss > (0.8 * ncell(x)))) {
				if (!is.null(ext)) {
					x <- crop(x, ext)
				}
				sr <- cbind(1:ncell(x), values(x))
				colnames(sr) <- c("cell", names(x))
			} else {
				sr <- spatSample(x, ss, "random", replace=replace, na.rm=TRUE, ext=ext, cells=TRUE, values=TRUE, warn=FALSE)
			}
		} else {
			if (!inherits(weights, "SpatRaster")) {
				error("spatSample", "weights must be a SpatRaster")
			}
			if (!compareGeom(x, weights)) {
				error("spatSample", "geometry of weights does not match the geometry of x")
			}
			sr <- vector("list", length=length(f))
			for (i in 1:length(f)) {
				r <- x == f[i]
				r <- mask(weights, r, maskvalue=TRUE, inverse=TRUE)
				sr[[i]] <- sampleWeights(r, size, replace=replace, values=FALSE, cells=TRUE, ext=ext)
			}
			sr <- unlist(sr)
			sr <- cbind(cell=sr, extract(x, sr))
		}
		uv <- unique(sr[,2])
		nuv <- length(uv)
		if (each) {
			sz <- rep(size, nuv)
		} else {
			sz <- rep(floor(size / nuv), nuv)
			d <- size - sum(sz)
			i <- sample(nuv, d)
			sz[i] <- sz[i] + 1
		}
		ys <- vector(mode="list", length=length(uv))
		for (i in seq_len(length(uv))) {
			y <- sr[sr[, 2] == uv[i], ,drop=FALSE]
			if (nrow(y) > sz[i]) {
				y <- y[sample.int(nrow(y), sz[i]),  ,drop=FALSE]
			}
			ys[[i]] <- y
		}
		res <- do.call(rbind, ys)
		colnames(res) <- c('cell', names(x))

		if (warn) {
			ta <- table(res[,2])
			sz <- sz[sz > 0]
			ta <- names(ta[ta < sz])
			tb <- f[!(f %in% unique(res[,2]))]
			tba <- c(tb, ta)
			if ((length(tba) > 0)) {
				warn("spatSample", "fewer samples than requested for group(s): ", paste(tba, collapse=", "))
			}
		}
	}

	if ((!xy) && (!as.points)) cells <- TRUE
	if (xy) {
		pts <- xyFromCell(x, res[,1])
		res <- cbind(res[,1,drop=FALSE], pts, res[,2,drop=FALSE])
	}
	if (!values) {
		res <- res[,1:(1 + 2*(xy|as.points)), drop=FALSE]
	}

	if (as.points) {
		if (!xy) {
			pts <- xyFromCell(x, res[,1])
		}
		res <- vect(pts, crs=crs(x), atts=data.frame(res))
	} else if (as.df) {
		res <- data.frame(res)
	}
	if (!cells) {
		res <- res[,-1,drop=FALSE]
	}
	res
}


.seed <- function() {
  sample.int(.Machine$integer.max, 1)
}


.sampleCellsMemory <- function(x, size, replace, lonlat, ext=NULL) {
	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	if (nlyr(x) > 1) {
		x <- subst(anyNA(x), 1, NA)
	}

	if (lonlat) {
		v <- cbind(cell=1:ncell(x), abs(cos(pi * values(init(x, "y")) / 360)), values(x))
		v <- v[!is.na(v[,3]),]
		ssize <- ifelse(replace, size, min(size, nrow(v)))
		i <- sample.int(nrow(v), ssize, prob=v[,2], replace=replace)
	} else {
		v <- cbind(cell=1:ncell(x), values(x))
		v <- v[!is.na(v[,2]),]
		ssize <- ifelse(replace, size, min(size, nrow(v)))
		i <- sample.int(nrow(v), ssize, replace=replace)
	}
	v[i,1]
}


.sampleCellsRandom <- function(x, size, replace, na.rm=FALSE, ext=NULL, exp=5) {
	r <- rast(x)
	lonlat <- is.lonlat(r, perhaps=TRUE, warn=TRUE)
	if (!is.null(ext)) {
		r <- crop(rast(r), ext)
	}

	if ((!replace) && (size >= ncell(r))) {
		cells <- 1:ncell(r)
	} else {
		if (na.rm) {
			esize <- size * exp
		} else {
			esize <- size
		}
		if (na.rm && (blocks(x, n=4)$n == 1)) {
			cells <- .sampleCellsMemory(x, esize, replace, lonlat, ext)
		} else if (lonlat) {
			m <- ifelse(replace, 1.5, 2)
			n <- m * esize
			y <- yFromRow(r, 1:nrow(r))
			rows <- sample.int(nrow(r), n, replace=TRUE, prob=abs(cos(pi*y/180)))
			cols <- sample.int(ncol(r), n, replace=TRUE)
			cells <- cellFromRowCol(r, rows, cols)
			if (!replace) {
				cells <- unique(cells)
			}
		} else {
			if (!replace) esize <- min(ncell(r), esize)
			cells <- sample.int(ncell(r), esize, replace=replace)
		}
	}
	if (!is.null(ext)) {
		cells <- cellFromXY(x, xyFromCell(r, cells))
	}
	if (na.rm) {
		v <- rowSums(is.na(x[cells])) == 0
		cells <- cells[v]
	}
	if (length(cells) > size) {
		cells <- cells[1:size]
	}
	cells
}

.sampleCellsStratified <- function(x, size, each=TRUE, replace=FALSE, ext=NULL) {
	if (!is.null(ext)) {
		r <- rast(x)
		if (window(x)) {
			ext <- intersect(ext(x), ext)
		}
		window(x) <- ext
	}
	s <- x@pntr$sampleStratifiedCells(size, each, replace, .seed(), spatOptions())
	s[[1]] <- s[[1]] + 1
	if (!is.null(ext)) {
		s[[1]] <- cellFromXY(x, xyFromCell(r, s[[1]]))
	}
	out <- do.call(cbind, s)
	colnames(out) <- c("cell", names(x)[1])
	out
}




.sampleCellsRegular <- function(x, size, ext=NULL, exact=FALSE) {
	r <- rast(x)
	lonlat <- is.lonlat(r, perhaps=TRUE, warn=TRUE)
	if (!is.null(ext)) {
		r <- crop(rast(r), ext)
	}

	if (prod(size) >= ncell(r)) {
		cells <- 1:ncell(r)
	} else {
		if (TRUE) { # skipping lonlat adjustment for now
			if (exact) {
				xy <- regular_exact(r, size)	
			} else {
				xy <- spatSample(ext(r), size, "regular", lonlat, FALSE)
			}
			cells <- cellFromXY(r, xy)		
		} else {
			if (lonlat) {
				#ratio <- ncol(r)/nrow(r)
				e <- ext(r)
				r1 = e$xmax - e$xmin;
				r2 = e$ymax - e$ymin;
				halfy = e$ymin + r2/2;

				# beware that -180 is the same as 180; and that latitude can only go from -90:90 therefore:
				dx = distance(cbind(e$xmin, halfy), cbind(e$xmin + 1, halfy), TRUE, TRUE) * min(180.0, r1);
				dy = distance(cbind(0, e$ymin), cbind(0, e$ymax), TRUE, TRUE);
				ratio = dy/dx;
				n <- sqrt(size)
				#nx <- max(1, (round(n*ratio)))
				#ny <- max(1, (round(n/ratio)))
				nx <- min(max(1, round(n/ratio)), ncol(r))
				ny <- min(max(1, round(n*ratio)), nrow(r))

				xi <- ncol(r) / nx
				yi <- nrow(r) / ny
				rows <- unique(round(seq(.5*yi, nrow(r), yi)))
				w <- cos(pi*yFromRow(r, rows)/180)
				w <- w * length(w)/sum(w)
				xi <- xi / w
				xi <- pmax(1,pmin(xi, ncol(r)))
				z <- list()

				# needs refinement:
				global <- diff(e[1:2]) > 355

				if (global) {
					xi <- round(ncol(r) / round(ncol(r) / xi))
					for (i in 1:length(rows)) {
						if (xi[i] == 1) {
							cols <- 1:ncol(r)
						} else {
							cols <- seq(xi[i]/2, ncol(r)-1, xi[i])
						}
						z[[i]] <- cbind(rows[i], cols)
					}
				} else {
	#				xi <- round(ncol(r) / (round((ncol(r) / xi))+1))
					for (i in 1:length(rows)) {
						cols <- seq(xi[i]/2, ncol(r), xi[i])
						z[[i]] <- cbind(rows[i], cols)
					}
				}
				z <- do.call(rbind, z)
				cells <- cellFromRowCol(r, z[,1], z[,2])
			} else {
				f <- sqrt(size / ncell(r))
				nr <- ceiling(nrow(r) * f)
				nc <- ceiling(ncol(r) * f);
				xstep <- ncol(r) / nc
				ystep <- nrow(r) / nr
				xsamp <- seq(0.5*xstep, ncol(r), xstep)
				ysamp <- seq(0.5*ystep, nrow(r), ystep)
				xy <- expand.grid(ysamp, xsamp)
				cells <- cellFromRowCol(r, xy[,1], xy[,2])
			}
		}
	}
	if (!is.null(ext)) {
		cells <- cellFromXY(x, xyFromCell(r, cells))
	}
	cells
}



set_factors <- function(x, ff, cts, asdf) {
	if (!asdf) return(x)
	
	if (any(ff)) {
		x <- data.frame(x)
		for (i in which(ff)) {
			ct <- cts[[i]]
			m <- match(x[[i]], ct[,1])
			if (!inherits(ct[[2]], "numeric")) {
				x[[i]] <- factor(ct[m,2], levels=unique(ct[[2]]))
			} else {
				x[[i]] <- ct[m,2]
			}
		}
	}
	data.frame(x)
}



.sampleCellsExhaustive <- function(x, size, replace, ext=NULL, weights=NULL, warn=TRUE) {

	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	rx <- rast(x)
	x <- cells(x)
	if (length(x) < size) {
		if (!replace) {
			warn("spatSample", "fewer samples than requested are available")
			return(x)
		}
		size <- length(x)
	}
	
	if (!is.null(weights)) {
		if (!inherits(weights, "SpatRaster")) {
			error("spatSample", "weights must be a SpatRaster")
		}
		weights <- weights[[1]]
		if (!is.null(ext)) {
			weights <- crop(weights, ext)
		}
		if (!compareGeom(x, weights)) {
			error("spatSample", "geometry of weights does not match the geometry of x")
		}
		weights <- weights[x]
		s <- sample.int(x, size, prob=weights, replace=replace)

	} else if (is.lonlat(rx)) {	
		y <- xyFromCell(rx, x)[,2]
		weights <- abs(cos(pi * y / 360))
		s <- sample(x, size, prob=weights, replace=replace)
	} else {
		s <- sample(x, size, replace=replace)
	}
	s
}

sampleRaster <- function(x, size, method, replace, ext=NULL, warn, overview=FALSE) {
#	hadWin <- hasWin <- FALSE
	if (!is.null(ext)) {
#		hasWin <- TRUE
#		hadWin <- window(x)
#		oldWin <- ext(x)
		w <- intersect(ext(x), ext(ext))
		if (is.null(w)) {
			error("sampleRaster", "x and ext do not intersect")
		}
		x <- deepcopy(x)
		window(x) <- w
	}
	if (method == "regular") {
		if (length(size) > 1) {
			x@pntr <- x@pntr$sampleRowColRaster(size[1], size[2], warn[1])
		} else {
			x@pntr <- x@pntr$sampleRegularRaster(size, overview)
		}
	} else if (method == "random") {
		x@pntr <- x@pntr$sampleRandomRaster(size, replace, .seed())
	} else {
		error("spatSample", "method must be 'regular' or 'random' if as.raster=TRUE")
	}
#	if (hasWin) {
#		window(x) <- NULL
#		if (any(hadWin)) {
#			window(x) <- oldWin
#		}
#	}
	messages(x, "spatSample")
}


add_cxyp <- function(x, cnrs, cells, xy, as.points, values, na.rm) {

	if (na.rm) {
		if (values) {
			vals <- x[cnrs]
			i <- rowSums(is.na(vals)) == 0
			vals <- vals[i, , drop=FALSE]
		} else if (hasValues(x)) {
			tmpvals <- x[cnrs]
			i <- rowSums(is.na(tmpvals)) == 0
		} 
		cnrs <- cnrs[i]
	} else if (values) {
		vals <- x[cnrs]
	}
	out <- NULL
	if (cells) {
		out <- matrix(cnrs, ncol=1)
		colnames(out) <- "cell"
	}
	if (xy) {
		out <- cbind(out, xyFromCell(x, cnrs))
	}
	if (values) {
		if (is.null(out)) {
			out <- vals
		} else {
			out <- cbind(out, vals)		
		}
	}
	if (as.points) {
		if (xy) {
			out <- data.frame(out)
			v <- vect(out, geom=c("x", "y"), crs=crs(x), keepgeom=TRUE)
		} else {
			crds <- xyFromCell(x, cnrs)
			# crds is a matrix, no geom argument
			v <- vect(crds, crs=crs(x))
			values(v) <- out
		}
		return(v)
	}
	out
}



sampleRandom <- function(x, size, replace=FALSE, na.rm=FALSE, as.raster=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL, warn=TRUE, weights=NULL, exp=5, exhaustive=FALSE) {

	ff <- is.factor(x)
	lv <- levels(x)
	size <- size[1]

	cvals <- NULL
	if (cells || xy || as.points) {
		if (exhaustive && na.rm) {
			cnrs <- .sampleCellsExhaustive(x, size, replace, ext, weights=NULL, warn=FALSE)
		} else {
			cnrs <- .sampleCellsRandom(x, size, replace, na.rm, ext, exp=exp)
		}
		
		if ((length(cnrs) < size) && warn) {
			warn("spatSample", "fewer cells returned than requested")
		} else if (length(cnrs) > size) {
			cnrs <- cnrs[1:size]
		}

		out <- add_cxyp(x, cnrs, cells, xy, as.points, values, na.rm)
		return(out)
	}
	if (!hasValues(x)) {
		error("spatSample", "SpatRaster has no values")
	}

	if (!is.null(ext)) {
		if (window(x)) {
			x <- crop(x, ext)
		} else {
			window(x) <- ext
		}
	}

	if (exhaustive && na.rm) {
		cnrs <- .sampleCellsExhaustive(x, size, replace, ext=NULL, weights=NULL, warn=FALSE)
		out <- x[cnrs]	
	} else {
		#v <- x@pntr$sampleRandomValues(size, replace, seed)
		
		if (size > 0.75 * ncell(x)) {
			if (na.rm) {
				out <- stats::na.omit(values(x))
				attr(out, "na.action") <- NULL
				if (nrow(out) < size) {
					if (replace) {
						out <- out[sample.int(nrow(out), size, replace=TRUE), ,drop=FALSE]
					} else {
						warn("spatSample", "more non-NA cells requested than available")
					}
				} else {
					out <- out[sample.int(nrow(out), size), ,drop=FALSE]
				}
			} else {
				out <- values(x)
				out <- out[sample.int(nrow(out), size, replace=replace), ,drop=FALSE]
			}
			out <- set_factors(out, ff, lv, as.df)
			return(out)
		} # else {
		if (na.rm) {
			scells <- NULL
			ssize <- size*2
			for (i in 1:10) {
				scells <- c(scells, .sampleCellsRandom(x, ssize, replace, na.rm))
				if ((i>1) && (!replace)) {
					scells <- unique(scells)
				}
				out <- extractCells(x, scells, raw=!as.df)   
				out <- stats::na.omit(out)
				if (nrow(out) >= size) {
					out <- out[1:size, ,drop=FALSE]
					attr(out, "na.action") <- NULL
					rownames(out) <- NULL
					break
				}
			}
		} else {
			scells <- .sampleCellsRandom(x, size, replace, na.rm=FALSE)
			out <- extractCells(x, scells, raw=!as.df)   
		}
		if (NROW(out) < size) {
			if (warn) warn("spatSample", "fewer values returned than requested")
		} else if (is.null(dim(out))) {
			out = out[1:size]
		} else {
			out = out[1:size, ,drop=FALSE]
		}
		#out <- set_factors(out, ff, lv, as.df)
		return(out)
	}
}


sampleStratified2 <- function(x, size, replace=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL, warn=TRUE, each=TRUE) {


	ff <- is.factor(x)
	size <- size[1]
	vc <- .sampleCellsStratified(x, size, each, replace, ext)
	

	if (each) {
		uc <- unique(vc[,2])
		if ((length(uc) * size) > nrow(vc)) {
			warn("spatSample", paste("not all classes have", size, "cells"))
		}
	} else if (size > nrow(vc)) {
		warn("spatSample", paste("not all classes had sufficient cells"))
	}

	cvals <- NULL
	if (cells || xy || as.points) {
		out <- add_cxyp(x, vc[,1], cells, xy, as.points, values=FALSE, na.rm=FALSE)
		if (values) out <- cbind(out, vc[,2,drop=FALSE])
		return(out)
	} else if (values) {
		if (any(ff) && as.df) {
			lv <- levels(x)
			vc[,2] <- set_factors(vc[,2], ff, lv, as.df)
		}
	} else {
		vc <- vc[,1,drop=FALSE]
	}
	return(vc)
}


sampleRegular <- function(x, size, replace=FALSE, na.rm=FALSE, as.raster=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL, exact=FALSE) {

	ff <- is.factor(x)
	lv <- levels(x)

	if (length(size) > 2) {
		error("spatSample", "size argument should have length 1 or 2")
	}

	if (length(size)==2) {
		exact <- FALSE
	}

	if (cells || xy || as.points) {
		if (exact) {
			pxy <- regular_exact(x, size)
			cnrs <- cellFromXY(x, pxy)
		} else if (length(size) == 2) {
			cnrs <- x@pntr$sampleRowCol(size[1], size[2])
		} else {
			cnrs <- .sampleCellsRegular(x, size, ext, exact=exact)
		}
		return(add_cxyp(x, cnrs, cells, xy, as.points, values, na.rm))
	}
	
	if (!hasValues(x)) {
		error("spatSample", "SpatRaster has no values")
	}

	if (!is.null(ext)) x <- crop(x, ext)

	if (exact && (length(size) == 1)) {
		xy <- regular_exact(x, size)
		v <- extract(x, xy, ID=FALSE)
	} else {
		opt <- spatOptions()
		if (length(size) > 1) {
			v <- x@pntr$sampleRowColValues(size[1], size[2], opt)
		} else {
			v <- x@pntr$sampleRegularValues(size, opt)
		}
		x <- messages(x, "spatSample")
			
		if (length(v) > 0) {
			v <- do.call(cbind, v)
			colnames(v) <- names(x)
		}
		v <- set_factors(v, ff, lv, as.df)
	}
	if (na.rm) {
		v <- na.omit(v)
	}
	return(v)
}


sampleRegularKM <- function(x, size, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL, ...) {

	ff <- is.factor(x)
	lv <- levels(x)
	cnrs <- .km_regular(x, size, ...)

	if (cells || xy || as.points) {
		return(add_cxyp(x, cnrs, cells, xy, as.points, values, na.rm=FALSE))
	}
	
	if (!hasValues(x)) {
		error("spatSample", "SpatRaster has no values")
	}

	if (!is.null(ext)) x <- crop(x, ext)
	x[cnrs]
}



setMethod("spatSample", signature(x="SpatRaster"),
	function(x, size, method="random", replace=FALSE, na.rm=FALSE, as.raster=FALSE, as.df=TRUE, as.points=FALSE, values=hasValues(x), cells=FALSE, xy=FALSE, ext=NULL, warn=TRUE, weights=NULL, exp=5, exhaustive=FALSE, exact=FALSE, each=TRUE, ...) {

		if (method == "display") return(sampleRaster(x, size, "regular", FALSE, ext=ext, warn=FALSE, overview=TRUE))
		method <- match.arg(tolower(method), c("random", "regular", "spread", "stratified", "weights"))

		if (!as.points) {
			if (!(values || cells || xy)) {
				error("spatSample", "at least one of 'values', 'cells', or 'xy' must be TRUE; or 'as.points' must be TRUE")
			}
		}

		exp <- max(c(1, exp), na.rm=TRUE)
		size <- round(size)
		if (isTRUE(any(size < 1)) || isTRUE(any(is.na(size)))) {
			error("spatSample", "sample size must be a positive integer")
		}

		if ((!replace) && (method != "regular")) {
			if (length(size) > 1) {
				error("spatSample", "sample size must be a single number")
			}
			if (warn && (size > ncell(x))) {
				warn("spatSample", "requested sample size is larger than the number of cells")
				size <- ncell(x)
			}
		}

		if (as.raster) return(sampleRaster(x, size, method, replace, ext, warn))

		lonlat <- is.lonlat(x, perhaps=TRUE, warn=FALSE)
		if (lonlat) exact <- FALSE

		if (method == "regular") {
			sampleRegular(x, size, replace=replace, na.rm=na.rm, as.df=as.df, as.points=as.points, values=values, cells=cells, xy=xy, ext=ext, exact=exact)
		} else if (method == "spread") {
			sampleRegularKM(x, size, as.df=as.df, as.points=as.points, values=values, cells=cells, xy=xy, ext=ext, ...)
		} else if (method == "stratified") {
			if (!is.null(weights)) {  # use old method
				return( sampleStratified_old(x, size, replace=replace, as.df=as.df, as.points=as.points, cells=cells, values=values, xy=xy, ext=ext, warn=warn, exp=exp, weights=weights, exhaustive=exhaustive, lonlat=lonlat, each=each) )
			} else {  # new method
				return( sampleStratified2(x, size, replace=replace, as.df=as.df, as.points=as.points, values=values, cells=cells, xy=xy, ext=ext, warn=warn, each=each) )			
			}
		} else if (!is.null(weights)) {  # should also implement for random
			error("spatSample", "argument weights is only used when method='stratified'")
		} else if (method == "random") {
			sampleRandom(x, size, replace=replace, na.rm=na.rm, as.df=as.df, as.points=as.points, values=values, cells=cells, xy=xy, ext=ext, warn=warn, exp=exp, exhaustive=exhaustive)
		} else if (method == "weights") {
			out <- try(sampleWeights(x, size, replace=replace, as.df=as.df, values=values, as.points=as.points, cells=cells, xy=xy, ext=ext) )
			if (inherits(out, "try-error")) {
				error("spatSample", "weighted sample failed. Perhaps the data set is too big")
			}
			return (out)
		}
	}
)


setMethod("spatSample", signature(x="SpatExtent"),
	function(x, size, method="random", lonlat, as.points=FALSE, exact=FALSE) {
		if (missing(lonlat)) {
			error("spatSample", "provide a lonlat argument")
		}
		if (lonlat) {
			stopifnot(x$ymax <= 90 || x$ymin >= -90)
		}
		method <- match.arg(method, c("regular", "random"))
		size <- round(size)
		stopifnot(size > 0)
		if (method=="random") {
			s <- x@pntr$sampleRandom(size, lonlat, .seed())
		} else if (exact) {
			s <- regular_exact(x, size)
			colnames(s) <- c("x", "y")
			if (as.points) {
				s <- vect(s)
			}
			return(s)
		} else {
			s <- x@pntr$sampleRegular(size, lonlat)
		}
		s <- do.call(cbind, s)
		colnames(s) <- c("x", "y")
		if (as.points) {
			s <- vect(s)
		}
		s
	}
)





.grid_sample <- function(xy, n=1, r, chess="") {

	cell <- cellFromXY(r, xy)
    uc <- unique(stats::na.omit(cell))

	chess <- trimws(chess)
	if (chess != "") {
		chess <- match.arg(tolower(chess), c("white", "black"))
		nc <- ncol(r)
		if (nc %% 2 == 1) {
			if (chess=="white") {
				tf <- 1:ceiling(ncell(r)/2) * 2 - 1
			} else {
				tf <- 1:ceiling((ncell(r)-1)/2) * 2
			}
		} else {
			nr <- nrow(r)
			row1 <- 1:(ceiling(nr / 2)) * 2 - 1
			row2 <- row1 + 1
			row2 <- row2[row2 <= nr]

			if (chess=="white") {
				col1 <- 1:(ceiling(nc / 2)) * 2 - 1
				col2 <- col1 + 1
				col2 <- col2[col2 <= nc]
			} else {
				col1 <- 1:(ceiling(nc / 2)) * 2
				col2 <- col1 - 1
				col1 <- col1[col1 <= nc]
			}

			cells1 <- cellFromRowColCombine(r, row1, col1)
			cells2 <- cellFromRowColCombine(r, row2, col2)
			tf <- c(cells1, cells2)
		}
		uc <- uc[uc %in% tf]
	}

    cell <- cellFromXY(r, xy)
    cell <- cbind(1:nrow(xy), cell, stats::runif(nrow(xy)))
	cell <- stats::na.omit(cell)

    cell <- cell[order(cell[,3]), ]
    sel <- list()
    for (i in 1:length(uc)) {
        ss <- subset(cell, cell[,2] == uc[i])
        sel[[i]] <- ss[1:min(n, nrow(ss)), 1]
    }
	unlist(sel)
}


#coordinates <- function(x) {
#	do.call(cbind, x@pntr$coordinates())
#}

get_field_name <- function(x, nms, sender="") {
	x <- x[1]
	if (is.numeric(x)) {
		x <- round(x)
		if (x > 0 && x <= length(nms)) {
			x = nms[x]
		} else {
			error(sender, "invalid index. there are ", length(nms), " columns")
		}
	} else if (is.character(x)) {
		if (!(x %in% nms)) {
			error(sender, "invalid name")
		}
	}
	x
}


setMethod("spatSample", signature(x="SpatVector"),
	function(x, size, method="random", strata=NULL, chess="") {
		method = match.arg(tolower(method), c("regular", "random"))
		size <- round(size)
		stopifnot(size > 0)
		gtype <- geomtype(x)
		if (gtype %in% c("lines", "polygons")) {
			if (!is.null(strata)) {
				if (length(strata) == 1) {
					if (is.character(strata)) {
						stopifnot(strata %in% names(x))
					} else  {
						stopifnot((strata > 0) && (strata < ncol(x)))
					}
					strata <- x[[strata, drop=TRUE]]
				} else if (length(strata) != length(x)) {
					stop("length of strata must be 1 or length(x)")
				}
				s <- stats::na.omit(unique(strata))
				n <- length(size)
				if (n==1) {
					n <- rep_len(n, length(s))
				} else if (length(s) != n) {
					stop("length of strata must be 1 or length(na.omit(unique(strata)))")
				}
				r <- lapply(s, function(s) {
					spatSample(x[strata == s, ], size, method, NULL, "")
				})
				r <- do.call(rbind, r)
				return(r)
			}
			out <- vect()
			if (length(size) == 1) {
				out@pntr <- x@pntr$sample(size, method[1], .seed())
			} else {
				out@pntr <- x@pntr$sampleGeom(size, method[1], .seed())	
			}
			messages(x, "spatSample")
			return(messages(out, "spatSample"))
		} else if (grepl(gtype, "points")) {
			if (!is.null(strata)) {
				if (inherits(strata, "SpatRaster")) {
					xy <- crds(x)
					i <- .grid_sample(xy, size[1], rast(strata), chess)
					return(x[i,])
				} else {
					error("spatSample", "not yet implemented for these strata")
				}
			} else {
				error("spatSample", "use `sample` to sample (point) geometries")
			}
		}
	}
)

