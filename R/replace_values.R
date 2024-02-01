# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3

setMethod("set.values", signature(x="SpatRaster"),
	function(x, cells, values, layer=0)  {

		#if (any(is.na(cells))) {
		#	warn("set.values", "cells should not be NA")
		#}

		if (is.character(layer)) {
			layer <- match(layer, names(x))
			if (any(is.na(layer))) {
				error("set.values", "invalid layer")
			}
		}
		layer <- round(layer)

		if (all(layer > 0)) {
			if (missing(cells) && missing(values)) {
				return(invisible(TRUE));
			}
			if (any(is.na(layer))) { error("set.values", "layers cannot be NA")}
			if (inherits(layer, "character")) {
				layer <- match(layer, names(x))
				if (any(is.na(layer))) { error("set.values", "invalid layer names")}
			}
			if (any((layer < 1) | (layer > nlyr(x))))  { error("set.values", "invalid layer numbers") }
			n <- length(layer)
			if (n > length(unique(layer)))  { error("set.values", "duplicated layers") }

			bylyr <- FALSE
			if (!is.null(dim(values))) {
				if (ncol(values) != n) {
					error("set.values", "ncol(values) does not match `length(layer)`")
				}
				bylyr <- TRUE
				#if (inherits(values, "data.frame")) {
				#	values <- as.matrix(values)
				#}
				values <- as.vector(values)
			}
			ok <- x@ptr$replaceCellValuesLayer(layer-1, cells-1, values, bylyr, spatOptions())
			messages(x)
			invisible(TRUE)
		} else {
			if (any(layer > 0)) {
				error("set.values", "some (but not all) layer numbers are < 1")
			}
			if (missing(cells) && missing(values)) {
				return(readAll(x))
			}
			bylyr <- FALSE
			if (!is.null(dim(values))) {
				if (ncol(values) != nlyr(x)) {
					error("set.values", "ncol(values) does not match the nlyr(x)")
				}
				bylyr <- TRUE
				#if (inherits(values, "data.frame")) {
				#	values <- as.matrix(values)
				#}
				values <- as.vector(values)
			}
			ok <- x@ptr$replaceCellValues(cells-1, values, bylyr, spatOptions())
			messages(x)
		}
		invisible(TRUE)
	}
)


make_replace_index <- function(v, vmx, nreps, name="i") {
	caller <- paste0("`[<-`(", name, ")")

	if (inherits(v, "SpatRaster")) {
		error(caller, paste("index", name, "cannot be a SpatRaster"))
	} 
	if (inherits(v, "SpatVector")) {
		error(caller, paste("index", name, "cannot be a SpatVector"))
	}
	if (inherits(v, "SpatExtent")) {
		error(caller, paste("index", name, "cannot be a SpatExtent"))
	}

	if (!is.numeric(v)) {

		if (NCOL(v) > 1) {
			error(caller, paste("index", name, "has multiple columns"))
		}
		if (inherits(v, "data.frame")) {
			v <- v[,1,drop=TRUE]
		} else if (inherits(v, "matrix")) {
			v <- as.vector(v)
		} 
		if (!is.vector(v)) {
			error(caller, paste("the type of index", name, "is unexpected:", class(v)[1]))
		}
		if (is.factor(v) || is.character(v)) {
			error(caller, paste("the type of index", name, "cannot be a factor or character"))
		} 
		if (is.logical(v)) {
			if (length(v) > vmx) {
				error(caller, paste("index", name, "is too long"))
			} 
			if (length(v) <= vmx) {
				v <- which(rep_len(v, vmx))
			}
		} else {
			v <- as.numeric(v)
		}
	}

	if (inherits(v, "matrix")) {
		if (ncol(v) == 1) {
			v <- v[,1]
		} else if (nrow(v) == 1) {
			v <- v[1,]
		} else {
			error(caller, paste("index", name, "has unexpected dimensions:", paste(dim(v), collapse=", ")))
		}
	}

	if (any(is.na(v))) {
		if (nreps > 1) {
			error(caller, "NAs are not allowed in subscripted assignments")
		} else {
			v <- v[!is.na(v)]
		}
	}
	#vv <- stats::na.omit(v)
	if (all(v < 0)) {
		if (any(v < -vmx)) {
			error(caller, paste(name, "is out of its valid range"))
		}
		v <- (1:vmx)[v]
	} 
	if (any(v < 1 | v > vmx)) {
		error(caller, paste(name, "is out of its valid range"))
	}
	v
}


.replace_all <- function(x, value) {
	nl <- nlyr(x)
	if (is.matrix(value)) {
		d <- dim(value)
		if (!all(d == c(ncell(x), nl))) {
			if ((d[2] == nl) && (d[1] < ncell(x))) {
				value <- apply(value, 2, function(i) rep_len(i, ncell(x)))
			} else {
				error("`[`","dimensions of the matrix do not match the SpatRaster")
			}
		}
		x <- try( setValues(x, value, TRUE, TRUE) )
	} else {
		v <- try( matrix(nrow=ncell(x), ncol=nl) )
		if (! inherits(x, "try-error")) {
			v[] <- value
			x <- try( setValues(x, v, TRUE, TRUE) )
		}
	}
	if (inherits(x, "try-error")) {
		error("`[`", "cannot set values")
	}
	return(x)
}

.replace_cell <- function(x, i, k, value) {
	bylyr = FALSE
	if (!is.null(dim(value))) {
		stopifnot(ncol(value) == nlyr(x))
		bylyr <- TRUE
		if (inherits(value, "data.frame")) {
			value <- as.matrix(value)
		}
		value <- as.vector(value)
	}
	opt <- spatOptions()
	x@ptr <- x@ptr$deepcopy()
	if (is.na(k[1])) {
		if (!x@ptr$replaceCellValues(i-1, value, bylyr, opt)) {
			messages(x, "`[<-`")
		} 
	} else {
		if (!x@ptr$replaceCellValuesLayer(k-1, i-1, value, bylyr, opt)) {
			messages(x, "`[<-`")
		}
	}
	x
}

.replace_cell_lyr <- function(x, cell, lyrs, value) {
	ulyrs <- sort(unique(lyrs))
	opt <- spatOptions()
	for (lyr in ulyrs) {
		y <- x[[lyr]]
		i <- which(lyrs == lyr)
		if (!y@ptr$replaceCellValues(cell[i]-1, value[i], FALSE, opt)) {
			messages(y, "`[<-`")
		}
		x[[lyr]] <- y
	}
	x
}


.replace_spatvector <- function(x, i, value) {
	if (length(value) > 1) {
		if (length(value) > nrow(i)) {
			# could be by layer if NCOL>1?
			error("`[`", "value is too long")
		}
		value <- rep_len(value, length.out=length(i))
	}
	rasterize(i, x, field=value, update=TRUE)
}


.replace_spatextent <- function(x, i, value) {
	if (length(value) > 1) {
		if (length(value) > nrow(i)) {
			# could be by layer if NCOL>1?
			error("`[`", "value is too long")
		}
		value <- rep_len(value, length.out=length(i))
	}
	rasterize(as.polygons(i), x, field=value, update=TRUE)
}


.replace_spatraster <- function(x, i, value) {
	if (inherits(value, "SpatRaster")) {
		x <- mask(x, i, maskvalues=TRUE)
		cover(x, value)
	} else {
		if (NCOL(value) > 1) {
			error("`[<-`", "cannot use a data.frame with multiple columns")
		}
		value <- unlist(value)
		if (length(value) == 1) {
			mask(x, i, maskvalues=TRUE, updatevalue=value[1])
		} else {
			i <- as.logical(values(i))
			i[is.na(i)] <- FALSE #TRUE, for #1115
			i <- which(i)
			x[i] <- value
			x
		}
	}
}



setReplaceMethod("[", c("SpatRaster", "ANY", "ANY", "ANY"),
	function(x, i, j, k, value) {

		m <- c(missing(i), missing(j), missing(k))
		s <- rep(FALSE, 3)
		if (!m[1]) s[1] <- inherits(i, "list")
		if (!m[2]) s[2] <- inherits(j, "list")
		if (!m[3]) s[3] <- inherits(k, "list")
		if (any(s)) {
			if (m[1]) i <- NULL 
			if (m[2]) j <- NULL 
			if (m[3]) k <- NULL 
			i <- rcl(x, i, j, k)
			m <- c(FALSE, TRUE, TRUE)
		}

		if (missing(value)) {
			value <- k
			k <- NA
			m[3] <- TRUE
		}

		if ((!m[1]) && (inherits(i, "matrix"))) {
			if (ncol(i) == 1) {
				i <- i[,1]
			} else if (ncol(i) == 2) {
				i <- cellFromRowCol(x, i[,1], i[,2])
				m[2] <- TRUE
				m[3] <- TRUE
			} else if (ncol(i) == 3) {
				k <- i[,3]
				value <- rep_len(value, length(k))
				i <- cellFromRowCol(x, i[,1], i[,2])
				return(.replace_cell_lyr(x, i, k, value))
			} else {
				error("`[<-`", paste("index i has", ncol(i), "columns"))
			}
		} 

		if (!m[3]) {
			if (inherits(k, "character")) {
				k <- match(k, names(x))
				if (any(is.na(k))) {
					stop()
				}
			} else {
				k <- make_replace_index(k, nlyr(x), length(value), "k")
			}
		} else {
			k <- NA
		}

		if (all(m)) {
			return(.replace_all(x, value))
		} 

		if (!m[1]) { # i not missing
			if (inherits(i, "SpatRaster")) {
				return(.replace_spatraster(x, i, value))
			} 
			if (inherits(i, "SpatVector")) {
				return(.replace_spatvector(x, i, value))
			}
			if (inherits(i, "SpatExtent")) {
				return(.replace_spatextent(x, i, value))
			}
			theCall <- sys.call(-1)
			narg <- length(theCall)-length(match.call(call=theCall))
			if ((narg==0) && m[2]) {
				# cell
				i <- make_replace_index(i, ncell(x), length(value), "i")
			} else if (m[2]) {
				# row
				i <- make_replace_index(i, nrow(x), length(value), "i")
				i <- cellFromRowColCombine(x, i, 1:ncol(x))
			} else {
				#row,col
				i <- make_replace_index(i, nrow(x), length(value), "i")
				j <- make_replace_index(j, ncol(x), length(value), "j")
				i <- cellFromRowColCombine(x, i, j)
			}
		} else if (!m[2]) {
			#col
			j <- make_replace_index(j, ncol(x), length(value), "j")
			i <- cellFromRowColCombine(x, 1:nrow(x), j)
		} else {
			if (inherits(value, "SpatRaster")) {
				x[[k]] <- value
				return(x)
			}
			i <- 1:ncell(x)
		}
		return(.replace_cell(x, i, k, value))
	} 
)

