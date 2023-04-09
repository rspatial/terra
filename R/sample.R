
sampleWeights <- function(x, size, replace=FALSE, as.df=TRUE, as.points=FALSE, cells=FALSE, xy=FALSE, ext=NULL) {
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
	res <- res[i,]
	if (as.points) {
		res <- vect(res, c("x", "y"), crs=crs(x))
		if (!xy) {
			res$x <- NULL
			res$y <- NULL
		}
	} else if (as.df) {
		res <- data.frame(res)
	}
	res
}


sampleStratMemory <- function(x, size, replace, lonlat, ext=NULL, weights=NULL, warn=TRUE) {
	if (!is.null(ext)) {
		x <- crop(x, ext)
		if (!is.null(weights)) {
			weights <- crop(weights, ext)
		}
	}
	
	if (!is.null(weights)) {
		if (!inherits(weights, "SpatRaster")) {
			error("spatSample", "weights must be a SpatRaster")
		}
		if (!compareGeom(x, weights)) {
			error("spatSample", "geometry of weights does not match the geometry of x")
		}
		v <- cbind(cell=1:ncell(x), values(x), values(weights))
		v <- v[!is.na(v[,2]), ]
		uv <- sort(unique(v[,2]))
		ys <- vector(mode="list", length=length(uv))
		for (i in seq_len(length(uv))) {
			vv <- v[v[,2] == uv[i], ]
			if (replace) {
				s <- sample.int(nrow(vv), size, prob=vv[,3], replace=TRUE)
			} else {
				s <- sample.int(nrow(vv), min(size, nrow(vv)), prob=vv[,3], replace=FALSE)
			}
			ys[[i]] <- vv[s,-3]
		}
	} else if (lonlat) {
		v <- cbind(cell=1:ncell(x), values(x), abs(cos(pi * values(init(x, "y")) / 360)))
		v <- v[!is.na(v[,2]), ]
		uv <- sort(unique(v[,2]))
		ys <- vector(mode="list", length=length(uv))
		for (i in seq_len(length(uv))) {
			vv <- v[v[,2] == uv[i], ]
			if (replace) {
				s <- sample.int(nrow(vv), size, prob=vv[,3], replace=TRUE)
			} else {
				s <- sample.int(nrow(vv), min(size, nrow(vv)), prob=vv[,3], replace=FALSE)
			}
			ys[[i]] <- vv[s,-3]
		}
	} else {
		v <- cbind(cell=1:ncell(x), values(x))
		v <- v[!is.na(v[,2]), ]
		uv <- sort(unique(v[,2]))
		ys <- vector(mode="list", length=length(uv))
		for (i in seq_len(length(uv))) {
			vv <- v[v[,2] == uv[i], ]
			if (replace) {
				s <- sample.int(nrow(vv), size, replace=TRUE)
			} else {
				s <- sample.int(nrow(vv), min(size, nrow(vv)), replace=FALSE)
			}
			ys[[i]] <- vv[s,]
		}
	}
	ys <- do.call(rbind, ys)
	if (warn) {
		ta <- tapply(ys[,1], ys[,2], length)
		ta <- names(ta)[ta < size]
		if (length(ta) > 0) {
			warn("spatSample", 'fewer samples than requested are available for group(s): ', paste(ta, collapse=', '))
		}
	}
	ys
}



sampleStratified <- function(x, size, replace=FALSE, as.df=TRUE, as.points=FALSE, cells=TRUE, xy=FALSE, ext=NULL, warn=TRUE, exp=5, weights=NULL) {

	if (nlyr(x) > 1) {
		x <- x[[1]]
		warn("spatSample", "only the first layer of x is used")
	}
	if (!hasValues(x)) {
		error("spatSample", "x has no values")
	}

	lonlat <- is.lonlat(x, perhaps=TRUE, warn=FALSE)
	if (blocks(x, n=4)$n == 1) {
		res <- sampleStratMemory(x, size, replace, lonlat, ext, weights, warn)
	} else {

		#f <- freq(x)
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
				sr[[i]] <- sampleWeights(r, size, replace=replace, cells=TRUE, ext=ext)[,1]
			}
			sr <- unlist(sr)
			sr <- cbind(cell=sr, extract(x, sr))
		}
		uv <- unique(sr[,2])
		ys <- vector(mode="list", length=length(uv))
		for (i in seq_len(length(uv))) {
			y <- sr[sr[, 2] == uv[i], ,drop=FALSE]
			if (nrow(y) > size) {
				y <- y[sample.int(nrow(y), size),  ,drop=FALSE]
			}
			ys[[i]] <- y
		}
		res <- do.call(rbind, ys)
		colnames(res) <- c('cell', names(x))

		if (warn) {
			ta <- table(res[,2])
			ta <- names(ta[ta < size])
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
		i <- sample.int(nrow(v), min(size, nrow(v)), prob=v[,2], replace=replace)
	} else {
		v <- cbind(cell=1:ncell(x), values(x))
		v <- v[!is.na(v[,2]),]
		i <- sample.int(nrow(v), min(size, nrow(v)), replace=replace)
	}
	v[i,1]
}


.sampleCells <- function(x, size, method, replace, na.rm=FALSE, ext=NULL, exp=5) {
	r <- rast(x)
	lonlat <- is.lonlat(r, perhaps=TRUE, warn=TRUE)
	if (!is.null(ext)) {
		r <- crop(rast(r), ext)
	}

	if ((!replace) && (size >= ncell(r))) {
		cells <- 1:ncell(r)
	} else if (method == "random") {
		if (na.rm) {
			esize <- size * exp
		} else {
			esize <- size
		}
		if (na.rm && (blocks(x, n=4)$n == 1)) {
			cells <- .sampleCellsMemory(x, size, replace, lonlat, ext)
		} else if (lonlat) {
			m <- ifelse(replace, 1.5, 1.25)
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
	} else { # regular
		if (TRUE) {
			xy <- spatSample(ext(r), size, method, lonlat, FALSE)
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
	if (na.rm) {
		v <- rowSums(is.na(x[cells])) == 0
		cells <- cells[v]
	}
	if (method == "random") {
		if (length(cells) > size) {
			cells <- cells[1:size]
		}
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

sampleRaster <- function(x, size, method, replace, ext=NULL, warn) {
#	hadWin <- hasWin <- FALSE
	if (!is.null(ext)) {
#		hasWin <- TRUE
#		hadWin <- window(x)
#		oldWin <- ext(x)
		w <- intersect(ext(x), ext(ext))		
		x <- deepcopy(x)
		window(x) <- w
	}
	if (method == "regular") {
		if (length(size) > 1) {
			x@pnt <- x@pnt$sampleRowColRaster(size[1], size[2], warn[1])
		} else {
			x@pnt <- x@pnt$sampleRegularRaster(size)
		}
	} else if (method == "random") {
		x@pnt <- x@pnt$sampleRandomRaster(size, replace, .seed())
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
	


setMethod("spatSample", signature(x="SpatRaster"),
	function(x, size, method="random", replace=FALSE, na.rm=FALSE, as.raster=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL, warn=TRUE, weights=NULL, exp=5, exhaustive=FALSE) {

		exp <- max(c(1, exp), na.rm=TRUE)
		size <- round(size)
		if (isTRUE(any(size < 1)) || isTRUE(any(is.na(size)))) {
			error("spatSample", "sample size must be a positive integer")
		}
		method <- match.arg(tolower(method), c("random", "regular", "stratified", "weights"))

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

		if (method == "stratified") {
			return( sampleStratified(x, size, replace=replace, as.df=as.df, as.points=as.points, cells=cells, xy=xy, ext=ext, warn=warn, exp=exp, weights=weights) )
		} else if (!is.null(weights)) {
			error("spatSample", "argument weights is only used when method='stratified'")
		}

		if (method == "weights") {
			if (nlyr(x) > 1) {
				x <- x[[1]]
				warn("spatSample", "only the first layer of x is used")
			}
			if (!hasValues(x)) {
				error("spatSample", "x has no values")
			}
			out <- try(sampleWeights(x, size, replace=replace, as.df=as.df, as.points=as.points, cells=cells, xy=xy, ext=ext) )
			if (inherits(out, "try-error")) {
				error("spatSample", "weighted sample failed. Perhaps the data set is too big")
			}
			return (out)
		}

		ff <- is.factor(x)
		lv <- levels(x)

		if (cells || xy || as.points) {

			size <- size[1]
			if (exhaustive && (method=="random") && na.rm) {
				cnrs <- .sampleCellsExhaustive(x, size, replace, ext, weights=NULL, warn=FALSE)
			} else {
				cnrs <- .sampleCells(x, size, method, replace, na.rm, ext)
			}
			
			if (method == "random") {
				if ((length(cnrs) < size) && warn) {
					warn("spatSample", "fewer cells returned than requested")
				} else if (length(cnrs) > size) {
					cnrs <- cnrs[1:size]
				}
			}
			out <- NULL
			
			if (cells) {
				out <- matrix(cnrs, ncol=1)
				colnames(out) <- "cell"
			}
			
			if (xy) {
				out <- cbind(out, xyFromCell(x, cnrs))
			}
			if (values && hasValues(x)) {
				e <- extract(x, cnrs)
				if (is.null(out)) {
					out <- e
				} else {
					out <- cbind(out, e)
				}
			}
			if (as.points) {
				if (xy) {
					out <- data.frame(out)
					out <- vect(out, geom=c("x", "y"), crs=crs(x))
				} else {
					xy <- xyFromCell(x, cnrs)
					# xy is a matrix, no geom argument
					v <- vect(xy, crs=crs(x))
					values(v) <- out
					return(v)
				}
			}
			return(out)
		}
		if (!hasValues(x)) {
			error("spatSample", "SpatRaster has no values")
		}

		#method <- tolower(method)
		#stopifnot(method %in% c("random", "regular"))
	
		if (!is.null(ext)) x <- crop(x, ext)

		if (method == "regular") {
			opt <- spatOptions()
			if (length(size) > 1) {
				v <- x@pnt$sampleRowColValues(size[1], size[2], opt)
			} else {
				v <- x@pnt$sampleRegularValues(size, opt)
			}
			x <- messages(x, "spatSample")
			if (length(v) > 0) {
				v <- do.call(cbind, v)
				colnames(v) <- names(x)
			}
			v <- set_factors(v, ff, lv, as.df)
			return(v)
		} else { # random
			size <- size[1]

			if (exhaustive && na.rm) {
				cnrs <- .sampleCellsExhaustive(x, size, replace, ext, weights=NULL, warn=FALSE)
				out <- x[cnrs]	
			} else {
				#v <- x@pnt$sampleRandomValues(size, replace, seed)
				if (size > 0.75 * ncell(x)) {
					if (na.rm) {
						out <- stats::na.omit(values(x))
						attr(x, "na.action") <- NULL
						if (nrow(out) < size) {
							warn("spatSample", "more non-NA cells requested than available")
						} else {
							out <- out[sample.int(nrow(out), size), ,drop=FALSE]
						}
					} else {
						out <- values(x)
						out <- out[sample.int(nrow(out), size, replace=replace), ,drop=FALSE]
					}
					out <- set_factors(out, ff, lv, as.df)
					return(out)
				}

				if (na.rm) {
					scells <- NULL
					ssize <- size*2
					for (i in 1:10) {
						scells <- c(scells, .sampleCells(x, ssize, method, replace, na.rm))
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
					scells <- .sampleCells(x, size, method, replace)
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
	}
)


setMethod("spatSample", signature(x="SpatExtent"),
	function(x, size, method="random", lonlat, as.points=FALSE) {
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
			s <- x@pnt$sampleRandom(size, lonlat, .seed())
		} else {
			s <- x@pnt$sampleRegular(size, lonlat)
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
#	do.call(cbind, x@pnt$coordinates())
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
		if (gtype == "polygons") {
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
				out@pnt <- x@pnt$sample(size, method[1], .seed())
			} else {
				out@pnt <- x@pnt$sampleGeom(size, method[1], .seed())	
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
		} else {
			error("spatSample", "not yet implemented for lines")
		}
	}
)

#spatSample(disagg(as.points(v)), 1, "stratified", strata=r, chess="")



# setMethod("spatSample", signature(x="SpatExtent"),
	# function(x, size, method="regular", lonlat, ...) {
		# if (missing(lonlat)) {
			# stop("provide a lonlat argument")
		# }
		# method = match.arg(method, c("regular", "random"))
		# size <- round(size)
		# stopifnot(size > 0)
		# e <- as.vector(x)
		# if (method=="random") {
			# if (lonlat) {
				# d <- round((e[4] - e[3]) * 1000);
				# dx <- (e[4] - e[3]) / (2 * d)
				# r <- unique(seq(e[3], e[4], length.out=d))
				# w <- abs(cos(pi*r/180))
				# x <- sample.int(length(r), size, prob=w, replace=TRUE)
				# lat <- r[x] + stats::runif(size, -dx, dx)
				# lon <- stats::runif(size, min = e[1], max = e[2])
				# vect(cbind(lon,lat), crs="+proj=lonlat +datum=WGS84")
			# } else {
				# x <- stats::runif(size, min = e[1], max = e[2])
				# y <- stats::runif(size, min = e[3], max = e[4])
				# vect(cbind(x, y))
			# }
		# } else {
			# r <- range(x)
			# ratio <- 0.5 * r[1]/r[2]
			# n <- sqrt(size)
			# nx <- max(1, (round(n*ratio)))
			# ny <- max(1, (round(n/ratio)))
			# xi <- r[1] / nx
			# yi <- r[2] / ny
			# if (lonlat) {
				# lat <- seq(e[3]+0.5*yi, e[4], yi)
				# w <- cos(pi*lat/180)
				# w <- w * length(w)/sum(w)
				# xi <- xi / w
				# xi <- pmin(xi, 180)
				# z <- list()
				# #off <- stats::runif(1)
				# for (i in 1:length(lat)) {
					# z[[i]] <- cbind(seq(e[1]+0.5*xi[i], e[2], xi[i]), lat[i])
				# }
				# z <- do.call(rbind, z)
				# vect(z, crs="+proj=lonlat +datum=WGS84")
			# } else {
				# x <- seq(e[1]+0.5*xi, e[2], xi)
				# y <- seq(e[3]+0.5*yi, e[4], yi)
				# vect(cbind(rep(x, length(y)), rep(y, each=length(x))))
			# }
		# }
	# }
# )
