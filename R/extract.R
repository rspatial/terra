# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# License GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warning("big number", "cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}


.makeDataFrame <- function(x, v) {
	v <- data.frame(v, check.names = FALSE)
	ff <- is.factor(x)
	if (any(ff)) {
		ff <- which(ff)
		levs <- levels(x)
		for (f in ff) {
			lev <- levs[[f]]
			v[[f]] = factor(v[[f]], levels=lev$levels)
			levels(v[[f]]) = lev$labels
		}
	}
	v
}

setlabs <- function(x, labs) {
	x[ (x<1) | (x>length(labs))] <- NA
	factor(labs[x], levels=labs)
}


wmean <- function(p) {
	n <- length(p)
	w <- p[[n]]
	p[[n]] <- NULL
	sapply(p, function(x) {
		stats::weighted.mean(x, w, na.rm = TRUE)
	})
}



setMethod("extract", signature(x="SpatRaster", y="SpatVector"), 
function(x, y, fun=NULL, method="simple", list=FALSE, factors=TRUE, cells=FALSE, xy=FALSE, weights=FALSE, touches=is.lines(y), ...) { 
	if (!is.null(fun)) {
		cells <- FALSE
		xy <- FALSE
	} 
	e <- x@ptr$extractVector(y@ptr, touches[1], method[1], isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]))
	x <- messages(x, "extract")


	#f <- function(i) if(length(i)==0) { NA } else { i }
	#e <- rapply(e, f, how="replace")
	cn <- names(x)
	if (!is.null(fun)) {
		if (weights) {
			test1 <- isTRUE(try( deparse(fun)[2] == 'UseMethod(\"mean\")', silent=TRUE))
			test2 <- isTRUE(try( fun@generic == "mean", silent=TRUE))
			if (!(test1 | test2)) { warn("extract", "the weighted mean is returned") }
			e <- t(sapply(e, wmean))
		} else {
			fun <- match.fun(fun) 
			e <- rapply(e, fun, ...)
		}
		e <- matrix(e, nrow=nrow(y), byrow=TRUE)
		colnames(e) <- cn
		e <- cbind(ID=1:nrow(e), e)
	} else {
		if (cells) {
			cn <- c(cn, "cell")
			i <- which(cn=="cell")
			e <- lapply(e, function(j) {
					j[[i]] <- j[[i]] + 1
					#names(j) <- cn
					j
				}
			)			
		}
		if (weights) cn <- c(cn, "weight")
		if (xy) {
			cn <- c(cn, "x", "y")
		}	
		if (!list) {
			e <- lapply(1:length(e), function(i) {
				ee <- unlist(e[[i]])
				if (length(ee) == 0) ee <- NA
				cbind(ID=i, matrix(ee, ncol=length(e[[i]])))
			})
			e <- do.call(rbind, e)
			colnames(e)[-1] <- cn	
		}
	}
	if (factors) {
		if (is.matrix(e)) {
			e <- data.frame(e, check.names = FALSE)
		}
		f <- is.factor(x)
		if (any(f)) {
			g <- cats(x)
			for (i in which(f)) {
				labs <- g[[i]]$labels
				if (!list) {
					v <- e[[i+1]] + 1
					if (!(is.null(fun))) {
						if (max(abs(v - trunc(v)), na.rm=T) > 0) {
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
	extract(x, y, ...)[,-1,drop=FALSE]
	#}
})


setMethod("extract", signature(x="SpatRaster", y="data.frame"), 
function(x, y, ...) { 
	if (ncol(y) != 2) {
		error("extract", "extract expects a 2 column data.frame of x and y coordinates")
	}
	extract(x, as.matrix(y), ...)
})


setMethod("extract", signature(x="SpatRaster", y="numeric"), 
function(x, y, ...) { 
	y <- as.integer(y)
	y[(y < 1) | (y > ncell(x))] <- NA
	x[y]
})


setMethod("[", c("SpatRaster", "missing", "missing"),
function(x, i, j, ... , drop=FALSE) {
	values(x, mat=!drop)
})

setMethod("[", c("SpatRaster", "logical", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x[which(i),, drop=drop]
})


extract_cell <- function(x, cells, drop=FALSE) {
	e <- x@ptr$extractCell(cells-1)
	messages(x, "[")
	e <- do.call(cbind, e)
	colnames(e) = names(x)
	.makeDataFrame(x, e)[,,drop]
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
	g <- geomtype(x)
	if (!grepl("points", g)) {
		stop("the first argument must be points")
	}
	r <- relate(x, y, "within")
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
		e <- cbind(1:nrow(x), e)
	}
	if (ncol(y) > 0) {
		d <- as.data.frame(y)	
		e <- data.frame(id.x=e[,1], d[e[,2], ,drop=FALSE])
		rownames(e) <- NULL
	} else {
		colnames(e) <- c("id.x", "id.y")
	}
	e
})

