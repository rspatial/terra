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
					error("set.values", "ncol(values) does not match the `length(layer)`")
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
				x@ptr$readAll()
				return(invisible(TRUE));
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


make_index <- function(v, vmx, name="i") {
	if (inherits(v, "SpatRaster")) {
		error("`[<-`", paste("index", name, "cannot be a SpatRaster"))
	} 
	if (inherits(v, "SpatVector")) {
		error("`[<-`", paste("index", name, "cannot be a SpatVector"))
	}

	if (!is.numeric(v)) {
		if (inherits(v, "data.frame")) {
			v <- v[,1,drop=TRUE]
		} 
		if (is.logical(v)) {
			if (length(v) > vmx) {
				error("`[<-`", paste(name, "is too long"))
			} 
			if (length(v) < vmx) {
				v <- which(rep_len(v, vmx))
			}
		} else {
			v <- as.numeric(v)
		}
	}
	if (any(v < 1 | v > vmx)) {
		error("`[<-`", paste(name, "is out of its valid range"))
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

.replace_cell <- function(x, i, value) {
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
	if (!x@ptr$replaceCellValues(i-1, value, bylyr, opt)) {
		messages(x, "`[<-`")
	} else {
		x
	}
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
			i[is.na(i)] <- TRUE
			i <- which(i)
			x[i] <- value
			x
		}
	}
}



setReplaceMethod("[", c("SpatRaster", "ANY", "ANY", "ANY"),
	function(x, i, j, k, value) {
		ni <- missing(i)
		nj <- missing(j)
		nk <- missing(k)
		if (missing(value)) {
			value = k
			k = NA
			nk <- TRUE
		} else if (!nk) {
			if (inherits(k, "character")) {
				k <- match(k, names(x))
				if (any(is.na(k))) {
					stop()
				}
			} else {
				k <- make_index(k, nlyr(x), "k")
			}
			if ((length(k) == nlyr(x)) && (all(k == 1:nlyr(x)))) {
				x <- .replace_all(x, value)
			} else {
			# should be able to pass "k" to the cpp methods instead
				y <- x[[k]]
				if (ni & nj) {
					y[] <- value 
				} else if (ni) {
					y[,j] <- value 
				} else {
					y[i,] <- value
				}
				x[[k]] <- y
			}
			return(x)
		}

		if (ni & nj) {
			return(.replace_all(x, value))
		} 
		
		if (!ni) { # i not missing			
			if (inherits(i, "SpatRaster")) {
				return(.replace_spatraster(x, i, value))
			} 
			if (inherits(i, "SpatVector")) {
				return(.replace_spatvector(x, i, value))
			}
			theCall <- sys.call(-1)
			narg <- length(theCall)-length(match.call(call=theCall))		
			if ((narg==0) && nj) {
				# cell
				i <- make_index(i, ncell(x), "i")
			} else if (nj) {
				# row
				i <- make_index(i, nrow(x), "i")
				i <- cellFromRowColCombine(x, i, 1:ncol(x))
			} else {
				#row,col
				i <- make_index(i, nrow(x), "i")
				j <- make_index(j, ncol(x), "j")
				i <- cellFromRowColCombine(x, i, j)
			}
		} else { #if (!nj) {
			#col
			j <- make_index(j, ncol(x), "j")
			i <- cellFromRowColCombine(x, 1:nrow(x), j)
		}
		return(.replace_cell(x, i, value))
	} 
)




