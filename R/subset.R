# Authors: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# License GPL v3

positive_indices <- function(i, n, na.rm=TRUE, caller="`[`") {
	i <- stats::na.omit(i)
	if (!is.logical(i)) {
		stopifnot(is.numeric(i))
		if (!(all(i <= 0) || all(i >= 0))) {
			error(caller, "you cannot mix positive and negative indices")
		}
	}
	i <- (1:n)[i]
	if (na.rm) {
		i[!is.na(i)]
	} else {
		i
	}
}


setMethod("subset", signature(x="SpatRaster"),
	function(x, subset, negate=FALSE, NSE=FALSE, filename="", overwrite=FALSE, ...) {
		if (NSE) {
			subset <- if (missing(subset)) {
				1:nlyr(x)
			} else {
				nl <- as.list(seq_along(names(x)))
				names(nl) <- nms <- names(x)
				v <- eval(substitute(subset), nl, parent.frame())
				if (!inherits(substitute(subset), "character")) {
					if (sum(nms %in% nms[v]) > length(v)) {
						error("subset", "you cannot select a layer with a name that is not unique")
					}
				}
				v
			}
		}
		if (is.character(subset)) {
			nms <- names(x)
			if (!all(subset %in% nms)) {
				error("subset", "invalid name(s)")
			}
			if (sum(nms %in% subset) > length(subset)) {
				error("subset", "you cannot select a layer with a name that is not unique")
			}
			subset <- match(subset, nms)
		}
		if (any(is.na(subset))) {
			error("subset", "undefined layer(s) selected")
		}
		if (negate) subset = -subset
		subset <- positive_indices(subset, nlyr(x), TRUE, "subset")
		opt <- spatOptions(filename, overwrite, ...)
		x@ptr <- x@ptr$subset(subset-1, opt)
		messages(x, "subset")
	}
)



## exact matching
setMethod("$", "SpatRaster",
	function(x, name) {
		subset(x, name, NSE=FALSE)
	}
)

setMethod("[[", c("SpatRaster", "character","missing"),
function(x, i, j) {
	subset(x, i, NSE=FALSE)
})


setMethod("[[", c("SpatRaster", "logical", "missing"),
function(x, i, j) {
	subset(x, which(i), NSE=FALSE)
})


setMethod("[[", c("SpatRaster", "numeric", "missing"),
function(x, i, j) {
	subset(x, i, NSE=FALSE)
})


setMethod("[[", c("SpatRaster", "ANY", "missing"),
function(x, i, j) {
	i <- as.vector(unlist(i))
	x[[i]]
})



setMethod("subset", signature(x="SpatVector"),
	function(x, subset, select, drop=FALSE, NSE=FALSE) {
 		if (NSE) {
			d <- as.list(x)
			# from the subset<data.frame> method
			r <- if (missing(subset)) {
					TRUE
				} else {
					r <- eval(substitute(subset), d, parent.frame())
					if (!is.logical(r)) error("subset", "argument 'subset' must be logical")
					r & !is.na(r)
				}
			v <- if (missing(select)) {
					TRUE
				} else {
					nl <- as.list(seq_along(d))
					names(nl) <- names(d)
					eval(substitute(select), nl, parent.frame())
				}
			x <- x[r, v, drop=drop]
		} else {
			if (missing(select)) {
				x <- x[which(as.vector(subset)), drop=drop]
			} else {
				x <- x[which(as.vector(subset)), select, drop=drop]
			}
		}
		#g <- gc()
		x
	}
)


.subset_cols <- function(x, subset, drop=FALSE) {
	if (is.character(subset)) {
		i <- stats::na.omit(match(subset, names(x)))
		if (length(i)==0) {
			error("subset", "no valid variable name found")
		} else if (length(i) < length(subset)) {
			warn("subset", "invalid variable name(s) excluded")
		}
	} else {
		i <- positive_indices(subset, ncol(x), TRUE, "subset")
		if (length(i)==0) {
			i <- 0
		}
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
	i <- positive_indices(i, nrow(x), TRUE, "`[`")
	x@ptr <- x@ptr$subset_rows(i-1)
	x <- messages(x, "[")
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})

setMethod("[", c("SpatVector", "numeric", "logical"),
function(x, i, j, drop=FALSE) {
	j <- which(rep_len(j, ncol(x)))
	x[i, j, drop=drop]
})


setMethod("[", c("SpatVector", "logical", "missing"),
function(x, i, j, drop=FALSE) {
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
function(x, i, j, drop=FALSE) {
	i <- positive_indices(i, nrow(x), TRUE, "`[`")
	j <- positive_indices(j, ncol(x), TRUE, "`[`")
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
function(x, i, j, drop=FALSE) {
	j <- positive_indices(j, ncol(x), TRUE, "`[`")
	x@ptr <- x@ptr$subset_cols(j-1)
	x <- messages(x, "[")
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})

setMethod("[", c("SpatVector", "missing", "character"),
function(x, i, j, drop=FALSE) {
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
function(x, i, j, drop=FALSE) {
	j <- which(rep_len(j, ncol(x)))
	x[,j,drop=drop]
})


setMethod("[", c("SpatVector", "numeric", "character"),
function(x, i, j, drop=FALSE) {
	j <- stats::na.omit(match(j, names(x)))
	if (length(j) == 0) j <- 0
	x <- x[i,j,drop=drop]
})

setMethod("[", c("SpatVector", "logical", "character"),
function(x, i, j, drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	x[i,j,drop=drop]
})


setMethod("[", c("SpatVector", "logical", "numeric"),
function(x, i, j, drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	x[i,j,drop=drop]
})

setMethod("[", c("SpatVector", "logical", "logical"),
function(x, i, j, drop=FALSE) {
	i <- which(rep_len(i, nrow(x)))
	j <- which(rep_len(j, ncol(x)))
	x[i,j,drop=drop]
})


setMethod("[", c("SpatVector", "missing", "missing"),
function(x, i, j,  drop=FALSE) {
	if (drop) {
		values(x)
	} else {
		x
	}
})


setMethod("[", c("SpatVector", "matrix", "missing"),
function(x, i, j, drop=FALSE) {
	if (ncol(i) > 1) warn("`[`", "using the first column of i to subset x")
	x[i[,1]]
})

setMethod("[", c("SpatVector", "data.frame", "missing"),
function(x, i, j, drop=FALSE) {
	if (ncol(i) > 1) warn("`[`", "using the first column of i to subset x")
	x[i[,1]]
})


setMethod("[", c("SpatVector", "data.frame", "ANY"),
function(x, i, j, drop=FALSE) {
	x <- x[i[,1], j]
})

setMethod("[", c("SpatVector", "character", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x[, i]
})


