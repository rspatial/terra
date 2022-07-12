# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# License GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warn("big number", "cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}


ext_from_rc <- function(x, r1, r2, c1, c2){
	e <- as.vector(ext(x))
	r <- res(x)
	c1 <- min(max(c1, 1), ncol(x))
	c2 <- min(max(c2, 1), ncol(x))
	if (c1 > c2) {
		tmp <- c1
		c1 <- c2
		c2 <- tmp
	}
	r1 <- min(max(r1, 1), nrow(x))
	r2 <- min(max(r2, 1), nrow(x))
	if (r1 > r2) {
		tmp <- r1
		r1 <- r2
		r2 <- tmp
	}

	xn <- xFromCol(x, c1) - 0.5 * r[1]
	xx <- xFromCol(x, c2) + 0.5 * r[1]
	yx <- yFromRow(x, r1) + 0.5 * r[2]
	yn <- yFromRow(x, r2) - 0.5 * r[2]
	ext(c(sort(c(xn, xx))), sort(c(yn, yx)))
}


#setlabs <- function(x, labs) {
#	x[ (x<1) | (x>length(labs))] <- NA
#	x <- factor(x, levels=1:length(labs))
#	levels(x) <- labs
#	x
#}


wmean <- function(p, na.rm=FALSE) {
	n <- length(p)
	w <- p[[n]]
	p[[n]] <- NULL
	sapply(p, function(x) {
		stats::weighted.mean(x, w, na.rm=na.rm)
	})
}

wsum <- function(p, na.rm=FALSE) {
	n <- length(p)
	w <- p[[n]]
	p[[n]] <- NULL
	sapply(p, function(x) {
		sum(x * w, na.rm=na.rm)
	})
}

wmin <- function(p, na.rm=FALSE) {
	n <- length(p)
	p[[n]] <- NULL
	sapply(p, function(x) {
		min(x, na.rm=na.rm)
	})
}

wmax <- function(p, na.rm=FALSE) {
	n <- length(p)
	p[[n]] <- NULL
	sapply(p, function(x) {
		max(x, na.rm=na.rm)
	})
}


		#if (!list) {
			#if (geomtype(y) == "points")  {
			#	e <- cbind(ID=1:length(e), matrix(unlist(e), ncol=nlyr(x), byrow=TRUE))
			#} else {
			#	e <- lapply(1:length(e), function(i) {
			#		ee <- unlist(e[[i]])
			#		if (length(ee) == 0) ee <- NA
			#		cbind(ID=i, matrix(ee, ncol=length(e[[i]])))
			#	})
			#	e <- do.call(rbind, e)
			#}
		#}



extractCells <- function(x, y, method="simple", list=FALSE, factors=TRUE, cells=FALSE, xy=FALSE, layer=NULL) {

	method <- match.arg(tolower(method), c("simple", "bilinear"))

	nl <- nlyr(x)
	useLyr <- FALSE
	if (!is.null(layer) && nl > 1) {
		if (any(is.na(layer))) {error("extract", "argument 'layer' cannot have NAs")}
		stopifnot(length(layer) == nrow(y))
		if (is.numeric(layer)) {
			layer <- round(layer)
			stopifnot(min(layer) > 0 & max(layer) <= nlyr(x))
		} else {
			layer <- match(layer, names(x))
			if (any(is.na(layer))) error("extract", "names in argument 'layer' do not match names(x)")
		}
		useLyr <- TRUE
	}
	cn <- names(x)
	opt <- spatOptions()
	if ((method == "bilinear") && (NCOL(y) > 1)) {
		e <- x@ptr$bilinearValues(y[,1], y[,2])
	} else {
		if (NCOL(y) == 2) {
			y <- cellFromXY(x, y)
		}
		e <- x@ptr$extractCell(y-1)
	}

	if (list) {
		messages(x, "extract")
		return(e)
	}
	e <- do.call(cbind, e)
	cn <- names(x)
	nc <- nl
	if (cells) {
		cn <- c(cn, "cell")
		nc <- nc + 1
		if (NCOL(y) == 2) {
			e <- cbind(e, cellFromXY(x, y))
		} else {
			e <- cbind(e, y)
		}
	}
	if (xy) {
		cn <- c(cn, "x", "y")
		nc <- nc + 2
		if (NCOL(y) == 1) {
			y <- xyFromCell(x, y)
		}
		e <- cbind(e, y)
	}
	colnames(e) <- cn

	if (factors) {
		if (method != "simple") {
			e <- as.data.frame(e)
		} else {
			e <- .makeDataFrame(x, e, TRUE)
		}
	}

	if (useLyr) {
		idx <- cbind(e[,1], layer[e[,1]]+1)
		ee <- cbind(e[,1,drop=FALSE], names(x)[idx[,2]-1], value=e[idx])
		colnames(ee)[2] <- "layer"
		if (ncol(e) > (nl+1)) {
			cbind(ee, e[,(nl+1):ncol(e), drop=FALSE])
		} else {
			ee
		}
	} else {
		e
	}
}

setMethod("extract", signature(x="SpatRaster", y="matrix"),
function(x, y, ...) {
	.checkXYnames(colnames(y))
	extractCells(x, y, ...)
})


setMethod("extract", signature(x="SpatRaster", y="SpatVector"),
function(x, y, fun=NULL, method="simple", list=FALSE, factors=TRUE, cells=FALSE, xy=FALSE, weights=FALSE, exact=FALSE, touches=is.lines(y), layer=NULL, ...) {

	nl <- nlyr(x)
	useLyr <- FALSE
	method <- match.arg(tolower(method), c("simple", "bilinear"))
	hasfun <- !is.null(fun)
	if (weights && exact) {
		exact = FALSE
	}
	if (hasfun) {
		cells <- FALSE
		xy <- FALSE
		if (weights || exact) {
			list <- TRUE
			fun <- .makeTextFun(fun)
			bad <- FALSE
			if (is.character(fun)) {
				if (!(fun %in% c("sum", "mean", "min", "max"))) {
					bad <- TRUE
				} else if (fun == "mean") {
					fun <- wmean
				} else if (fun == "sum") {
					fun <- wsum
				} else if (fun == "min") {
					fun <- wmin
				} else if (fun == "max") {
					fun <- wmax
				}
			} else {
				bad <- TRUE
			}
			if (bad) error("extract", 'if weights or exact=TRUE, "fun" must be "sum", "mean", "min", or "max"')
		}
	}
	if (!is.null(layer) && nl > 1) {
		if (any(is.na(layer))) {error("extract", "argument 'layer' cannot have NAs")}
		if (length(layer) == 1) {
			lyr_name <- layer
			layer <- as.character(y[[layer,drop=TRUE]])
		} else {
			lyr_name <- "layer"
			stopifnot(length(layer) == nrow(y))
		}
		if (is.numeric(layer)) {
			layer <- round(layer)
			stopifnot(min(layer) > 0 & max(layer) <= nl)
		} else {
			layer <- match(layer, names(x))
			if (any(is.na(layer))) error("extract", "names in argument 'layer' do not match names(x)")
		}
		useLyr <- TRUE
	}

	#f <- function(i) if(length(i)==0) { NA } else { i }
	#e <- rapply(e, f, how="replace")
	cn <- names(x)
	opt <- spatOptions()
	if (list) {
		e <- x@ptr$extractVector(y@ptr, touches[1], method, isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]), isTRUE(exact[1]), opt)
		x <- messages(x, "extract")
		if (weights || exact) {
			if (hasfun) {
				e <- sapply(e, fun, ...)
				e <- matrix(e, nrow=nrow(y), byrow=TRUE)
				colnames(e) <- cn
				e <- cbind(ID=1:nrow(e), e)
			}
		}
		return(e)
	}

	e <- x@ptr$extractVectorFlat(y@ptr, touches[1], method, isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]), isTRUE(exact[1]), opt)
	x <- messages(x, "extract")
	nc <- nl
	if (cells) {
		cn <- c(cn, "cell")
		nc <- nc + 1
	}
	if (weights) {
		cn <- c(cn, "weight")
		nc <- nc + 1
	} else if (exact) {
		cn <- c(cn, "fraction")
		nc <- nc + 1
	}
	if (xy) {
		cn <- c(cn, "x", "y")
		nc <- nc + 2
	}

	geo <- geomtype(y)
	if (geo == "points") {
		## this was? should be fixed upstream
		if (nc == nl) {
			e <- matrix(e, ncol=nc)
		} else {
			e <- matrix(e, ncol=nc, byrow=TRUE)
		}
		e <- cbind(1:nrow(e), e)
		if (nrow(e) > nrow(y)) { #multipoint
			g <- geom(y)
			e[,1] <- g[,1]
		}
	} else {
		e <- matrix(e, ncol=nc+1, byrow=TRUE)
	}
	cn <- c("ID", cn)
	colnames(e) <- cn
	if (hasfun) {
		fun <- match.fun(fun)
		e <- data.frame(e)
		e <- aggregate(e[,-1,drop=FALSE], e[,1,drop=FALSE], fun, ...)

		m <- sapply(e, NCOL)
		if (any(m > 1)) {
			e <- do.call(cbind, as.list(e))
			skip <- (length(cn) - nlyr(x))
			nms <- colnames(e)
			snms <- nms[(skip+1):length(nms)]
			mr <- max(m)
			if (!all(snms=="")) {
				snms <- paste0(rep(names(x), each=mr), ".", snms)
			} else {
				snms <- paste0(rep(names(x), each=mr), ".", rep(1:mr))
			}
			snms <- c(cn[1:skip], snms)
			colnames(e) <- snms
			e <- data.frame(e)
		}
	} else if (cells) {
		cncell <- cn =="cell"
		e[, cncell] <- e[, cncell] + 1
	}

	if (factors) {
		if (hasfun || method != "simple") {
			e <- as.data.frame(e)
		} else {
			id <- data.frame(e[,1,drop=FALSE])
			e <- cbind(id, .makeDataFrame(x, e[,-1,drop=FALSE], TRUE))
		}
	}

	if (useLyr) {
		idx <- cbind(e[,1], layer[e[,1]]+1)
		ee <- cbind(e[,1,drop=FALSE], names(x)[idx[,2]-1], value=e[idx])
		colnames(ee)[2] <- lyr_name
		if (ncol(e) > (nl+1)) {
			e <- cbind(ee, e[,(nl+1):ncol(e), drop=FALSE])
		} else {
			e <- ee
		}
	}
	e
})

setMethod("extract", signature(x="SpatRaster", y="sf"),
	function(x, y, fun=NULL, method="simple", list=FALSE, factors=TRUE, cells=FALSE, xy=FALSE, weights=FALSE, exact=FALSE, touches=is.lines(y), layer=NULL, ...) {
		y <- vect(y)
		extract(x, y, fun=fun, method=method, list=list, factors=factors, cells=cells, xy=xy, weights=weights, exact=exact, touches=touches, layer=layer, ...)
	}
)



setMethod("[", c("SpatRaster", "SpatVector", "missing"),
function(x, i, j, ... , drop=FALSE) {
	v <- extract(x, i)
	if (drop) {
		as.vector(v)
	} else {
		v
	}
})

setMethod("[", c("SpatVector", "SpatVector", "missing"),
function(x, i, j, ... , drop=FALSE) {
	r <- !relate(x, i, "disjoint")
	r <- which(apply(r, 1, any))
	x[r, ]
})


setMethod("[", c("SpatVector", "SpatExtent", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x[as.polygons(i)]
})


setMethod("extract", signature(x="SpatRaster", y="data.frame"),
function(x, y, ...) {
	if (ncol(y) != 2) {
		error("extract", "extract expects a 2 column data.frame of x and y coordinates")
	}
	v <- vect(y, colnames(y))
	extract(x, v, ...)
})


setMethod("extract", signature(x="SpatRaster", y="numeric"),
function(x, y, ...) {
	y <- round(y)
	y[(y < 1) | (y > ncell(x))] <- NA
	extract_cell(x, y)	
})

setMethod("extract", signature(x="SpatRaster", y="SpatExtent"),
function(x, y, factors=TRUE, cells=FALSE, xy=FALSE) {
	y <- cells(x, y)
	if (factors) dataframe = TRUE
	v <- extract_cell(x, y, factors=factors)
	if (cells) {
		v$cell <- y
	}
	if (xy) {
		v <- cbind(v, xyFromCell(x, y))
	}
	v
}
)


setMethod("[", c("SpatRaster", "missing", "missing"),
function(x, i, j, ... , drop=FALSE) {
	values(x, mat=!drop)
})

setMethod("[", c("SpatRaster", "logical", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x[which(i),, drop=drop]
})


extract_cell <- function(x, cells, drop=FALSE, factors=TRUE) {
	e <- x@ptr$extractCell(cells-1)
	messages(x, "extract_cell")
	e <- do.call(cbind, e)
	colnames(e) <- names(x)
	.makeDataFrame(x, e, factors)[,,drop]
}


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... , drop=TRUE) {

	add <- any(grepl("drop", names(match.call())))
	if (!drop) {
		if (nargs() == 3) {
			rc <- rowColFromCell(x, i)
			e <- ext_from_rc(x, min(rc[,1]), max(rc[,1]), min(rc[,2]), max(rc[,2]))
		} else {
			e <- ext_from_rc(x, min(i), max(i), 1, ncol(x))
		}
		return(crop(x, e))
	}
	if (nargs() > (2+add)) {
		i <- cellFromRowColCombine(x, i, 1:ncol(x))
	}
	extract_cell(x, i, drop=FALSE)
})


setMethod("[", c("SpatRaster", "data.frame", "missing"),
function(x, i, j, ... , drop=TRUE) {
	if (ncol(i) == 1) {
		i <- i[,1]
	} else if (ncol(i) == 2) {
		i <- cellFromXY(x, i)
	} else {
		error(" [", "cannot extract values with this data.frame")
	}
	`[`(x, i, drop=drop)
})

setMethod("[", c("SpatRaster", "matrix", "missing"),
function(x, i, j, ... , drop=TRUE) {
	if (ncol(i) == 1) {
		i <- i[,1]
	} else if ((nrow(i) == 1) && (ncol(i) != 2)) {
		i <- i[1,]
	} else if (ncol(i) == 2) {
		i <- cellFromXY(x, i)
	} else {
		error(" [", "cannot extract values with a ` of these dimensions")
	}
	`[`(x, i, drop=drop)
})


setMethod("[", c("SpatRaster", "missing", "numeric"),
function(x, i, j, ... , drop=TRUE) {
	if (!drop) {
		e <- ext_from_rc(x, 1, nrow(x), min(j), max(j))
		return(crop(x, e))
	}

	i <- cellFromRowColCombine(x, 1:nrow(x), j)
	extract_cell(x, i, drop=FALSE)
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ..., drop=TRUE) {
	if (!drop) {
		e <- ext_from_rc(x, min(i), max(i), min(j), max(j))
		return(crop(x, e))
	}
	i <- cellFromRowColCombine(x, i, j)
	extract_cell(x, i, drop=FALSE)
})


setMethod("[", c("SpatRaster", "SpatRaster", "missing"),
function(x, i, j, ..., drop=TRUE) {

	if (!compareGeom(x, i, crs=FALSE, stopOnError=FALSE)) {
		return (x[ext(i), drop=drop])
	}
	if (drop) {
		if (is.bool(i)) {
			i <- as.logical(values(i))
		} else {
			i <- !is.na(values(i))
		}
		values(x)[i,]
	} else {
		if (is.bool(i)) {
			mask(x, i, maskvalues=FALSE)
		} else {
			mask(x, i)
		}
	}
})


setMethod("[", c("SpatRaster", "SpatExtent", "missing"),
function(x, i, j, ..., drop=FALSE) {
	x <- crop(x, i)
	if (drop) {
		values(x)
	} else {
		x
	}
})


setMethod("extract", c("SpatVector", "SpatVector"),
function(x, y, ...) {

	#r <- relate(y, x, "within")
	#e <- apply(r, 1, which)
	r <- relate(x, y, "covers")
	e <- apply(r, 2, which)
	if (length(e) == 0) {
		e <- list(e)
	}
	if (is.list(e)) {
		e <- lapply(1:length(e), function(i) {
			if (length(e[[i]]) == 0) {
				cbind(i, NA)
			} else {
				cbind(i, e[[i]])
			}
		})
		e <- do.call(rbind, e)
	} else {
		e <- cbind(1:nrow(y), e)
	}
	if (ncol(x) > 0) {
		d <- as.data.frame(x)
		e <- data.frame(id.y=e[,1], d[e[,2], ,drop=FALSE])
		rownames(e) <- NULL
	} else {
		colnames(e) <- c("id.y", "id.x")
	}
	e
})

