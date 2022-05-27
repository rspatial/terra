# Authors: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# License GPL v3

positive_indices <- function(i, n, caller=" [ ") {
	if (!(all(i <= 0) || all(i >= 0))) {
		error(caller, "you cannot mix positive and negative indices")
	}
	i <- stats::na.omit(i)
	(1:n)[i]
}


setMethod("subset", signature(x="SpatRaster"), 
function(x, subset, filename="", overwrite=FALSE, ...) {
	if (is.character(subset)) {
		i <- match(subset, names(x))
	} else {
		i <- as.integer(subset)
		i[(i<1) | (i>nlyr(x))] <- NA
	}
	if (any(is.na(i))) {
		error("subset", paste("undefined layer(s) selected:", paste(subset[is.na(i)], collapse=", ")))
	}
	opt <- spatOptions(filename, overwrite, ...)
	x@ptr <- x@ptr$subset(i-1, opt)
	messages(x, "subset")
	return(x)
} )

setMethod("[", c("SpatRaster", "SpatVector", "missing"),
	function(x, i, j, ... ,drop=TRUE) {
		if (drop) {
			extract(x, i, data.frame=TRUE)[ , -1, drop=FALSE]
		} else {
			crop(x, i, mask=TRUE)
		}
	}
)


## expression matching
setMethod("[", c("SpatRaster", "character", "missing"),
	function(x, i, j, ... ,drop=TRUE) {
		i <- grep(i, names(x))
		subset(x, i, ...)
	}
)

## exact matching

setMethod("[[", c("SpatRaster", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})

setMethod("$", "SpatRaster",  
	function(x, name) { 
		subset(x, name) 
	} 
)

setMethod("[[", c("SpatRaster", "logical", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, which(i), ...)
})


setMethod("[[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	i <- positive_indices(i, nlyr(x), " [[ ")
	subset(x, i, ...)
})



setMethod("subset", signature(x="SpatVector"), 
	function(x, subset, select, drop=FALSE) {
 		d <- values(x)
		if (missing(subset)) {
			r <- rep_len(TRUE, nrow(d))
		} else {
			r <- eval(substitute(subset), d, parent.frame())
			if (!is.logical(r)) stop("'subset' must be logical")
			r <- r & !is.na(r)
		}
		if (missing(select)) 
			vars <- rep_len(TRUE, ncol(d))
		else {
			nl <- as.list(seq_along(d))
			names(nl) <- names(d)
			vars <- eval(substitute(select), nl, parent.frame())
		}
		x[r, vars, drop=drop]
	}
)


.subset_cols <- function(x, subset, drop=FALSE) {
	if (is.character(subset)) {
		i <- stats::na.omit(match(subset, names(x)))
	} else {
		i <- positive_indices(subset, ncol(x), "subset")
	}
	if (length(i)==0) {
		i <- 0
	} 
	if (length(i) < length(subset)) {
		warn(" [ ", "invalid columns omitted")
	}
	x@ptr <- x@ptr$subset_cols(i-1)
	x <- messages(x, "subset")
	if (drop) {	# drop geometry
		.getSpatDF(x@ptr$df)
	} else {
		x
	}
}


setMethod("[", c("SpatVector", "numeric", "missing"),
function(x, i, j, ... , drop=FALSE) {
	i <- positive_indices(i, nrow(x), "'['")
	x@ptr <- x@ptr$subset_rows(i-1)
	x <- messages(x, "[")
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})

setMethod("[", c("SpatVector", "numeric", "logical"),
function(x, i, j, ... , drop=FALSE) {
	j <- which(rep_len(j, ncol(x)))
	x[i, j, drop=drop]
})


setMethod("[", c("SpatVector", "logical", "missing"),
function(x, i, j, ... , drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	x@ptr <- x@ptr$subset_rows(i-1)
	x <- messages(x, "[")
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})

setMethod("[", c("SpatVector", "numeric", "numeric"),
function(x, i, j, ... , drop=FALSE) {
	i <- positive_indices(i, nrow(x), "'['")
	j <- positive_indices(j, ncol(x), "'['")
	p <- x@ptr$subset_rows(i-1)
	x@ptr <- p$subset_cols(j-1)
	x <- messages(x, "'['")
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})


setMethod("[", c("SpatVector", "missing", "numeric"),
function(x, i, j, ... , drop=FALSE) {
	j <- positive_indices(j, ncol(x), "'['")
	x@ptr <- x@ptr$subset_cols(j-1)
	x <- messages(x, "[")
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})

setMethod("[", c("SpatVector", "missing", "character"),
function(x, i, j, ... , drop=FALSE) {
	if (j[1] == "") {
		jj <- 0
	} else {
		jj <- match(j, names(x))
		if (any(is.na(jj))) {
			mis <- paste(j[is.na(jj)], collapse=", ")
			error(" x[,j] ", paste("name(s) not in x:", mis))
		}
		if (length(jj) == 0) { 
			jj <- 0
		}
	}
	x[,jj,drop=drop]
})

setMethod("[", c("SpatVector", "missing", "logical"),
function(x, i, j, ... , drop=FALSE) {
	j <- which(rep_len(j, ncol(x)))
	x[,j,drop=drop]
})


setMethod("[", c("SpatVector", "numeric", "character"),
function(x, i, j, ... , drop=FALSE) {
	j <- stats::na.omit(match(j, names(x)))
	if (length(j) == 0) j <- 0
	x <- x[i,j,drop=drop]
})

setMethod("[", c("SpatVector", "logical", "character"),
function(x, i, j, ... , drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	x[i,j,drop=drop]
})


setMethod("[", c("SpatVector", "logical", "numeric"),
function(x, i, j, ... , drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	x[i,j,drop=drop]
})

setMethod("[", c("SpatVector", "logical", "logical"),
function(x, i, j, ... , drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	j <- which(rep_len(j, ncol(x)))
	x[i,j,drop=drop]
})


setMethod("[", c("SpatVector", "missing", "missing"),
function(x, i, j, ... , drop=FALSE) {
	if (drop) {
		values(x)
	} else {
		x
	}
})


setMethod("[", c("SpatVector", "matrix", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x[i[,1]]
})

setMethod("[", c("SpatVector", "data.frame", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x[i[,1]]
})
