# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# License GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warn("big number", "cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}


.makeDataFrame <- function(x, v, factors=TRUE) {
	v <- data.frame(v, check.names = FALSE)
	if (factors) {
		ff <- is.factor(x)
		if (any(ff)) {
			ff <- which(ff)
			levs <- levels(x)
			for (f in ff) {
				lvs <- levs[[f]]
				v[[f]] = factor(v[[f]], levels=(1:length(lvs))-1)
				levels(v[[f]]) = levs[[f]]
			}
		}
	}
	v
}

setlabs <- function(x, labs) {
	x[ (x<1) | (x>length(labs))] <- NA
	x <- factor(x, levels=1:length(labs))
	levels(x) <- labs
	x
}


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
			if (bad) error("extract", 'if weights or exact=TRUE, "fun" must be "sum", "mean" or NULL')
		}
	} 
	if (!is.null(layer) && nl > 1) {
		if (any(is.na(layer))) {error("extract", "argument 'layer' cannot have NAs")}
		if (length(layer) == 1) {
			layer <- as.character(y[[layer,drop=TRUE]])
		} else {
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
	if (list) {
		e <- x@ptr$extractVector(y@ptr, touches[1], method, isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]), isTRUE(exact[1]))
		x <- messages(x, "extract")
		if (weights || exact) {
			if (hasfun) {
				e <- t(sapply(e, fun, ...))
				#e <- matrix(e, nrow=nrow(y), byrow=TRUE)
				colnames(e) <- cn
				e <- cbind(ID=1:nrow(e), e)
			}
		}
		return(e)
	}

	e <- x@ptr$extractVectorFlat(y@ptr, touches[1], method, isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]), isTRUE(exact[1]))
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
		if (is.matrix(e)) {
			e <- data.frame(e, check.names = FALSE)
		}
		f <- is.factor(x)
		if (any(f)) {
			g <- levels(x)
			for (i in which(f)) {
				labs <- g[[i]]
				if (!list) {
					v <- e[[i+1]] + 1
					if (!(is.null(fun))) {
						if (max(abs(v - trunc(v)), na.rm=TRUE) > 0) {
							next
						}
					}
					e[[i+1]] <- setlabs(v, labs)
				} else {
					e[[i]] <- rapply(e[[i]], function(j) setlabs(j+1, labs), how="replace")
				}
			}
		}
	}

	if (useLyr) {
		idx <- cbind(e[,1], layer[e[,1]]+1)
		ee <- cbind(e[,1,drop=FALSE], layer=names(x)[layer], value=e[idx])
		if (ncol(e) > (nl+1)) {
			e <- cbind(ee, e[,(nl+1):ncol(e), drop=FALSE])
		} else {
			e <- ee
		}
	}
	e
})


setMethod("[", c("SpatRaster", "SpatVector", "missing"),
function(x, i, j, ... , drop=FALSE) {
	v <- extract(x, i)
	if (drop) {
		as.vector(v)
	} else {
		v
	}
})


setMethod("extract", signature(x="SpatRaster", y="matrix"), 
function(x, y, ...) { 
	.checkXYnames(colnames(y))
	#if (length(list(...)) == 0) {
	#	i <- cellFromXY(x, y)
	#	r <- cbind(1:length(i), x[i])
	#	colnames(r) <- c("ID", names(x))
	#} else {
	y <- vect(y)
	e <- extract(x, y, ...)
	if (NCOL(e) > 1) {
		e[, -1, drop=FALSE]
	} else {
		e
	}
	#}
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
	y <- as.integer(y)
	y[(y < 1) | (y > ncell(x))] <- NA
	x[y]
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
	colnames(e) = names(x)
	.makeDataFrame(x, e, factors)[,,drop]
}

setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	if (nargs() > 2) {
		i <- cellFromRowColCombine(x, i, 1:ncol(x))
	} 
	extract_cell(x, i, drop)
})

setMethod("[", c("SpatRaster", "missing", "numeric"),
function(x, i, j, ... ,drop=FALSE) {
	i <- cellFromRowColCombine(x, 1:nrow(x), j)
	extract_cell(x, i, drop)
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ..., drop=FALSE) {
	i <- cellFromRowColCombine(x, i, j)
	extract_cell(x, i, drop)
})


setMethod("[", c("SpatRaster", "SpatRaster", "missing"),
function(x, i, j, ..., drop=FALSE) {
	x[which(as.logical(values(i))), drop=drop]
})


setMethod("extract", c("SpatVector", "SpatVector"),
function(x, y, ...) {
	r <- relate(y, x, "within")
	e <- apply(r, 1, which)
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
		e <- data.frame(id.x=e[,1], d[e[,2], ,drop=FALSE])
		rownames(e) <- NULL
	} else {
		colnames(e) <- c("id.x", "id.y")
	}
	e
})

