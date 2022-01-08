

sampleStratified <- function(x, size, replace=FALSE, as.df=TRUE, as.points=FALSE, cells=TRUE, xy=FALSE, ext=NULL, warn=TRUE, exp=2) {
	
	if ((!xy && !as.points)) cells <- TRUE
	
	f <- freq(x)	
	exp <- max(1, exp)
	ss <- exp * size * nrow(f)
	lonlat <- is.lonlat(x, perhaps=TRUE, warn=FALSE)
	if ((!lonlat) && (ss > (0.8 * ncell(x)))) { 
		sr <- cbind(1:ncell(x), values(x))
		colnames(sr) <- c("cell", names(x))
	} else {
		sr <- spatSample(x, ss, "random", replace=replace, na.rm=TRUE, ext=ext, cells=TRUE, values=TRUE, warn=warn)
	}
	
	ys <- list()
	notfound <- NULL

	for (i in seq_len(nrow(f))) {
		y <- sr[sr[, 2] == f[i,2], ,drop=FALSE]
		if (nrow(y) == 0) {
			notfound <- c(notfound, i)
		} else {
			if (nrow(y) > size) {
				y <- y[sample(nrow(y), size),  ,drop=FALSE]
			} 
			ys[[i]] <- y
		}
	}
	res <- do.call(rbind, ys)
	colnames(res) <- c('cell', names(x))
	
	ures <- unique(res[,2])
	miss <- !(ures %in% f[,"value"])
	if (any(miss) && warn) {
		miss <- which(miss)
		if (length(miss)== 1) {
			warn("sample", 'no samples for stratum: ', tanm)
		} else if (length(miss) > 1) {
			warn("sample", 'no samples for strata: ', paste(tanm, collapse=', '))
		}
	}
	
	ta <- tapply(res[,1], res[,2], length) 
	tanm <- names(ta)[which(ta < size)]
	if ((length(tanm) > 0) && warn) {
		if (length(tanm)== 1) {
			warn("sample", 'fewer samples than requested for stratum: ', tanm)
		} else if (length(tanm) > 1) {
			warn("sample", 'fewer samples than requested for strata: ', paste(tanm, collapse=', '))
		}
	}

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


.sampleCells <- function(x, size, method, replace, na.rm=FALSE, ext=NULL) {
	r <- rast(x)
	lonlat <- is.lonlat(r, perhaps=TRUE, warn=TRUE)
	if (!is.null(ext)) {
		r <- crop(rast(r), ext)
	}
	if ( ((!replace) || (method == "regular")) && (size >= ncell(r)) ) {
		cells <- 1:ncell(r)
	}

	if (method == "random") {
		nsize <- size
		if (na.rm) {
			if (replace) {
				size <- size*5
			} else {
				size <- min(ncell(r)*2, size*5)
			}
		}
		if (lonlat) {
			m <- ifelse(replace, 1.5, 1.25)
			n <- m * size
			y <- yFromRow(r, 1:nrow(r))
			w <- abs(cos(pi*y/180))
			rows <- sample.int(nrow(r), n, replace=TRUE, prob=w)
			cols <- sample.int(ncol(r), n, replace=TRUE)
			cells <- cellFromRowCol(r, rows, cols)
			if (!replace) {
				cells <- unique(cells)
			}
		} else {
			cells <- sample(ncell(r), size, replace=replace)
		}
	} else { # regular 
		if (lonlat) {
			ratio <- 0.5 * ncol(r)/nrow(r)
			n <- sqrt(size)
			nx <- max(1, (round(n*ratio)))
			ny <- max(1, (round(n/ratio)))
			xi <- ncol(r) / nx
			yi <- nrow(r) / ny
			rows <- unique(round(seq(.5*yi, nrow(r), yi)))

			w <- cos(pi*yFromRow(r, rows)/180)
			w <- w * length(w)/sum(w)
			xi <- xi / w
			xi <- pmax(1,pmin(xi, ncol(r)))
			z <- list()
			#off <- stats::runif(1) 
			for (i in 1:length(rows)) {
				z[[i]] <- cbind(rows[i], unique(round(seq(0.5*xi[i], ncol(r), xi[i]))))
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
	if (!is.null(ext)) {
		cells <- cellFromXY(x, xyFromCell(r, cells))
	}
	if (na.rm) {
		v <- rowSums(is.na(x[cells])) == 0
		cells <- cells[v]
	}
	if (method == "random") {
		if (length(cells) > nsize) {
			cells <- cells[1:nsize]
		}
	}
	return(cells)
}


set_factors <- function(x, ff, cts, asdf) {
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
	} else if (asdf) {
		x <- data.frame(x)
	}
	x
}


setMethod("spatSample", signature(x="SpatRaster"), 
	function(x, size, method="random", replace=FALSE, na.rm=FALSE, as.raster=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, xy=FALSE, ext=NULL, warn=TRUE) {

		size <- round(size)
		if (any(size < 1)) {
			error("spatSample", "sample size must be a positive integer")
		}
		if ((size > ncell(x)) & (!replace)) {
			size <- ncell(x)
		}

		method <- match.arg(tolower(method), c("random", "regular", "stratified"))
		if (method == "stratified") {
			if (as.raster) {
				error("as.raster is not valid for method='stratified'")
			}
			if (nlyr(x) > 1) {
				x <- x[[1]]
				warn("only the first layer of x is used")			
			}
			if (!hasValues(x)) {
				error("x has no values")			
			}
			return( sampleStratified(x, size, replace=replace, as.df=as.df, as.points=as.points, cells=cells, xy=xy, ext=ext, warn=warn, exp=5) )
		}

		if (!as.raster) {
			ff <- is.factor(x)
			lv <- active_cats(x)
		}

		if (cells || xy || as.points) {
			size <- size[1]
			cnrs <- .sampleCells(x, size, method, replace, na.rm, ext)
			if (method == "random") {
				if (length(cnrs) < size) {
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
				e <- set_factors(e, ff, lv, as.df)
				if (is.null(out)) {
					out <- e
				} else {
					out <- cbind(out, e)
				}
			}
			if (as.points) {
				if (xy) {
					out <- vect(out, crs=crs(x))
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
		if (!hasValues(x) & !as.raster) {
			error("spatSample", "SpatRaster has no values")
		}

		method <- tolower(method)
		stopifnot(method %in% c("random", "regular"))
		if (!replace) size <- pmin(ncell(x), size)

		if (!is.null(ext)) x <- crop(x, ext)

		if (method == "regular") {
			if (as.raster) {
				if (length(size) > 1) {
					x@ptr <- x@ptr$sampleRowColRaster(size[1], size[2])
				} else {
					x@ptr <- x@ptr$sampleRegularRaster(size)
				}
				x <- messages(x, "spatSample")
				return(x);
			} else {
				opt <- spatOptions()
				if (length(size) > 1) {
					v <- x@ptr$sampleRowColValues(size[1], size[2], opt)
				} else {
					v <- x@ptr$sampleRegularValues(size, opt)
				}
				x <- messages(x, "spatSample")
				if (length(v) > 0) {
					v <- do.call(cbind, v)
					colnames(v) <- names(x)
				}
				v <- set_factors(v, ff, lv, as.df)
				return(v)
			}
		} else { # random
			size <- size[1]
			if (as.raster) {
				x@ptr <- x@ptr$sampleRandomRaster(size, replace, .seed())
				x <- messages(x, "spatSample")
				return(x);
			} else {
				#v <- x@ptr$sampleRandomValues(size, replace, seed)
				if (size > 0.75 * ncell(x)) {
					if (na.rm) {
						out <- stats::na.omit(values(x))
						attr(x, "na.action") <- NULL
						if (nrow(out) < size) {
							warn("spatSample", "more non-NA cells requested than available")
						} else {
							out <- out[sample(nrow(out), size), ,drop=FALSE]
						}
					} else {
						out <- values(x)
						out <- out[sample(nrow(out), size, replace=replace), ,drop=FALSE]
					}
					out <- set_factors(out, ff, lv, as.df)
					return(out)
				}

				if (na.rm) {
					scells <- NULL
					ssize <- size*2
					for (i in 1:10) {
						scells <- c(scells, .sampleCells(x, ssize, method, replace))
						if ((i>1) && (!replace)) {
							scells <- unique(scells)
						}
						out <- stats::na.omit(x[scells])
						if (nrow(out) >= size) {
							out <- out[1:size, ,drop=FALSE]
							attr(out, "na.action") <- NULL
							rownames(out) <- NULL
							break
						}
					}
				} else {
					scells <- .sampleCells(x, size, method, replace)
					out <- x[scells]
				}
				if (NROW(out) < size) {
					if (warn) warn("spatSample", "fewer values returned than requested")
				} else if (method == "random") {
					if (is.null(dim(out))) {
						out = out[1:size]
					} else {
						out = out[1:size, ,drop=FALSE]
					}
				}
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
			s <- x@ptr$sampleRandom(size, lonlat, .seed())
		} else {
			s <- x@ptr$sampleRegular(size, lonlat)
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

	chess <- trim(chess)
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
#	do.call(cbind, x@ptr$coordinates())
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
			if (length(size) == 1) {
				x@ptr = x@ptr$sample(size, method[1], .seed())
			} else {
				x@ptr = x@ptr$sampleGeom(size, method[1], .seed())
			}
			return(messages(x))
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
