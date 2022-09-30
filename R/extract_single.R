
make_extract_index <- function(v, vmx, name="i") {

	caller <- paste0("`[`(", name, ")`")
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
		if (inherits(v, "data.frame")) {
			if (ncol(v) == 1) {
				v <- v[,1,drop=TRUE]
			} else if ((name == "i") && (ncol(v) == 2)) {
				v <- cellFromRowCol(x, v[,1], v[,2])
			} else {
				error(caller, paste("index", name, "has", ncol(v), "columns"))
			}
		} else if (inherits(v, "matrix")) {
			if (ncol(v) == 1) {
				v <- v[,1]
			} else {
				error(caller, paste("index", name, "is not numeric and has", ncol(v), "columns"))
			}
		}
		
		if (!is.vector(v)) {
			error(caller, paste("the type of index", name, "is unexpected:", class(v)[1]))		
		}
		if (is.factor(v) || is.character(v)) {
			error(caller, paste("the type of index", name, "cannot be a factor or character"))			
		}
		if (is.logical(v)) {
			v <- which(v)
		} 
		if (!is.numeric(v)) {
			error(caller, paste("the type of index", name, "is unexpected:", class(v)[1]))		
		}
	}
	if (inherits(v, "matrix")) {
		if (ncol(v) == 1) {
			v <- v[,1]
		} else {
			error(caller, paste("index", name, "has unexpected dimensions:", paste(dim(v), collapse=", ")))	
		}
	}
	positive_indices(v, vmx, caller=caller)
}



.extract_spatraster <- function(x, i, drop) {
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
}


.extract_spatextent <- function(x, i, drop) {
	x <- crop(x, i)
	if (drop) {
		values(x)
	} else {
		x
	}
}



.extract_spatvector <- function(x, i, drop) {
	if (drop) {
		extract(x, i, data.frame=TRUE)[ , -1, drop=FALSE]
	} else {
		crop(x, i, mask=TRUE)
	}
}



.extract_row <- function(x, i, drop=TRUE) {
	if (!drop) {
		e <- ext_from_rc(x, min(i), max(i), 1, ncol(x))
		return(crop(x, e))
	}
	i <- cellFromRowColCombine(x, i, 1:ncol(x))
	.extract_cell(x, i, drop=TRUE)
}


.extract_col <- function(x, j, drop=TRUE) {
	if (!drop) {
		e <- ext_from_rc(x, 1, nrow(x), min(j), max(j))
		return(crop(x, e))
	}
	i <- cellFromRowColCombine(x, 1:nrow(x), j)
	.extract_cell(x, i, drop=TRUE)
}

.extract_rowcol <- function(x, i, j, drop=TRUE) {
	if (!drop) {
		e <- ext_from_rc(x, min(i), max(i), min(j), max(j))
		return(crop(x, e))
	}
	i <- cellFromRowColCombine(x, i, j)
	.extract_cell(x, i, drop=TRUE)
}



.extract_cell <- function(x, i, drop=TRUE) {
	if (!drop) {	
		rc <- rowColFromCell(x, i)
		e <- ext_from_rc(x, min(rc[,1]), max(rc[,1]), min(rc[,2]), max(rc[,2]))
		crop(x, e)	
	} else {	
		e <- x@ptr$extractCell(i-1)
		messages(x, "extract")
		e <- do.call(cbind, e)
		colnames(e) <- names(x)
		.makeDataFrame(x, e)
	}
}	


.extract_cell_layer <- function(x, i, lyrs) {
	e <- x@ptr$extractCell(i-1)
	messages(x, "extract")
	e <- do.call(cbind, e)
	colnames(e) <- names(x)
	e <- .makeDataFrame(x, e)
	e[cbind(1:nrow(e), lyrs)]
}	



setMethod("[", c("SpatRaster", "ANY", "ANY", "ANY"),
	function(x, i, j, k, drop=TRUE) {
		if (missing(i)) {
			ni <- TRUE
			li <- FALSE
			i <- NULL
		} else {
			ni <- FALSE
			li <- is.list(i)
		}

		if (missing(j)) {
			nj <- TRUE
			lj <- FALSE
			j <- NULL
		} else {
			nj <- FALSE
			lj <- is.list(j)
		}

		if (missing(k)) {
			nk <- TRUE
			lk <- FALSE
			k <- NULL
		} else {
			nk <- FALSE
			lk <- is.list(k)
		}
		if (any(li, lj, lk)) {
			i <- rcl(x, i, j, k)
			nj <- nk <- TRUE
			ni <- FALSE
		}

			
		if (!nk) {
			if (is.logical(k) && length(k) == 1) {
				drop <- k
				nk <- TRUE
			}
		}
		
		if ((!ni) && (inherits(i, "matrix"))) {
			if (ncol(i) == 1) {
				i <- i[,1]
			} else if (ncol(i) == 2) {
				i <- cellFromRowCol(x, i[,1], i[,2])
				nj <- nk <- TRUE
			} else if (ncol(i) == 3) {
				k <- i[,3]
				i <- cellFromRowCol(x, i[,1], i[,2])
				return(.extract_cell_layer(x, i, k))
			} else {
				error("`[<-`", paste("index i has", ncol(i), "columns"))
			}
		} 		
		
		if (!nk) {
			if (is.logical(k) && length(k) == 1) {
				drop <- k
				nk <- TRUE
			} else {
				if (inherits(k, "character")) {
					k <- match(k, names(x))
					if (any(is.na(k))) {
						error("`[`(k)", "invalid layer name(s)")
					}
				} else {
					k <- make_extract_index(k, nlyr(x), "k")
				}
				x <- x[[k]]
			}
		}
		if ((!ni) && (inherits(i, "character"))) {
			# partial matching of layer names
			i <- grep(i, names(x))
			x <- subset(x, i, NSE=FALSE)
			if (nj) return(x)
			ni <- TRUE
		}


		if (ni && nj) {
			return(values(x, mat=!drop))
		} 
		
		if (!ni) { # i not missing			
			if (inherits(i, "SpatRaster")) {
				return(.extract_spatraster(x, i, drop))
			} 
			if (inherits(i, "SpatVector")) {
				return(.extract_spatextent(x, i, drop))
			}
			if (inherits(i, "SpatVector")) {
				return(.extract_spatextent(x, i, drop))
			}
			
			theCall <- sys.call(-1)
			narg <- length(theCall)-length(match.call(call=theCall))		
			if ((narg==0) && nj) {
				# cell
				i <- make_extract_index(i, ncell(x), "i")
				return(.extract_cell(x, i, drop=drop))
			} else if (nj) {
				# row
				i <- make_extract_index(i, nrow(x), "i")
				return(.extract_row(x, i, drop=drop))
			} else {
				#row,col
				i <- make_extract_index(i, nrow(x), "i")
				j <- make_extract_index(j, ncol(x), "j")
				return(.extract_rowcol(x, i, j, drop=drop))
			}
		} else { #if (!nj) {
			#col
			j <- make_extract_index(j, ncol(x), "j")
			return(.extract_col(x, j, drop=drop))
		}
	} 
)





setMethod("[", c("SpatVector", "SpatVector", "missing"),
function(x, i, j) {
	#r <- !relate(x, i, "disjoint")
	#r <- which(apply(r, 1, any))
	r <- is.related(x, i, "intersects")
	x[r, ]
})


setMethod("[", c("SpatVector", "SpatExtent", "missing"),
function(x, i, j) {
	x[as.polygons(i)]
})

