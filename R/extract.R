# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# License GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warning("big number", "cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}


.makeDataFrame <- function(x, v) {
	v <- data.frame(v)	
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


setMethod("extract", signature(x="SpatRaster", y="SpatVector"), 
function(x, y, fun=NULL, ..., touches=is.lines(y), method="simple", list=FALSE, factors=TRUE) { 
	e <- x@ptr$extractVector(y@ptr, touches[1], method[1])
	x <- messages(x, "extract")
	#f <- function(i) if(length(i)==0) { NA } else { i }
	#e <- rapply(e, f, how="replace")
	if (!is.null(fun)) {
		fun <- match.fun(fun) 
	  	e <- rapply(e, fun, ...)
		e <- matrix(e, nrow=nrow(y), byrow=TRUE)
		colnames(e) <- names(x)
		e <- cbind(ID=1:nrow(e), e)
	} else if (!list) {
		e <- lapply(1:length(e), function(i) cbind(ID=i, matrix(unlist(e[[i]]), ncol=length(e[[i]]))))
		e <- do.call(rbind, e)
		colnames(e)[-1] <- names(x)	
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
						if (max(abs(a - trunc(a)), na.rm=T) > 0) {
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
	if (length(list(...)) == 0) {
		i <- cellFromXY(x, y)
		r <- cbind(1:length(i), x[i])
		colnames(r) <- c("ID", names(x))
		r
	} else {
		y <- vect(y)
		extract(x, y, ...)
	}
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
	y[y < 1] <- NA
	y[y > ncell(x)] <- NA
	x[y]
})


setMethod("[", c("SpatRaster", "missing", "missing"),
function(x, i, j, ... , drop=FALSE) {
	values(x, mat=!drop)
})

setMethod("[", c("SpatRaster", "logical", "missing"),
function(x, i, j, ... , drop=FALSE) {
	v <- values(x)[as.logical(i), ]
	if (drop) {
		as.vector(v)
	} else {
		v
	}
})


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	if (any(stats::na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	if (nargs() > 2) {
		i <- cellFromRowCol(x, i, 1:ncol(x))
		# probably better to do return( readValues(x, i-1) )
	} 
	i[i<1] <- NA
	r <- x@ptr$extractCell(i-1)
	messages(x, "[")
	if (drop) {
		r
	} else {
		r <- do.call(cbind, r)
		colnames(r) = names(x)
		r
	}
})

setMethod("[", c("SpatRaster", "missing", "numeric"),
function(x, i, j, ... ,drop=FALSE) {
	i <- cellFromRowCol(x, 1:nrow(x), j)
	if (any(stats::na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	
	r <- x@ptr$extractCell(i-1)
	messages(x, "[")
	if (drop) {
		r
	} else {
		r <- do.call(cbind, r)
		colnames(r) = names(x)
		r
	}
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ..., drop=FALSE) {
	i <- cellFromRowColCombine(x, i, j)
	if (any(stats::na.omit(i) > 2^.Machine$double.digits)) .big_number_warning()
	r <- x@ptr$extractCell(i-1)
	messages(x, "[")
	if (drop) {
		r
	} else {
		r <- do.call(cbind, r)
		colnames(r) = names(x)
		r
	}
})



setMethod("[", c("SpatRaster", "SpatRaster", "missing"),
function(x, i, j, ..., drop=FALSE) {
	x[which(as.logical(values(i)))]
})

