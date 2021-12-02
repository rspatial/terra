

setMethod("$<-", "SpatRaster",  
	function(x, name, value) { 
		if (inherits(value, "SpatRaster")) {
			value <- value[[1]]
			names(value) <- name
		} else if (!is.null(value)) {
			y <- rast(x, nlyrs=1)
			test <- try(values(y) <- value, silent=TRUE)
			if (inherits(test, "try-error")) {
				error("$<-,SpatRaster", "the replacement value is not valid")
			}
			value <- y
			names(value) <- name
		}

		i <- which(name == names(x))[1]
		if (is.na(i)) {
			c(x, value)
		} else if (nlyr(x) == 1) {
			value
		} else if (i == 1) {
			c(value, x[[2:nlyr(x)]])
		} else if (i == nlyr(x)) {
			c(x[[1:(nlyr(x)-1)]], value)
		} else {
			c(x[[1:(i-1)]], value, x[[(i+1):nlyr(x)]])
		}
	}
)


setReplaceMethod("[[", c("SpatRaster", "character", "missing"),
	function(x, i, j, value) {
		if (inherits(value, "SpatRaster")) {
			if (nlyr(value) != length(i)) {
				error(" [[,SpatRaster", "length of names must be equal to the number of layers")
			}
			names(value) <- i
		} else if (length(i) > 1) {
			if (NCOL(value) > 1) {
				value <- as.list(data.frame(value))
			} else {
				stopifnot(length(i) == length(value))
			}
		} else if (!is.list(value)) {
			value <- list(value)
		}
		for (k in 1:length(i)) {
			eval(parse(text = paste0("x$", i[k], " <- value[[k]]")))
		}
		x
	}
)

setReplaceMethod("[[", c("SpatRaster", "numeric", "missing"),
	function(x, i, j, value) {
		if (nlyr(value) != length(i)) {
			error(" [[,SpatRaster,numeric", "length of indices must be equal to the number of layers")
		}
		if (any(i<1) | any(i > nlyr(x))) {
			error(" [[,SpatRaster,numeric", "indices must be between 1 and the number of layers")
		}
		for (k in 1:length(i)) {
			if (i[k] == 1) {
				x <- c(value[[k]], x[[2:nlyr(x)]])
			} else if (i[k] == nlyr(x)) {
				x <- c(x[[1:(nlyr(x)-1)]], value[[k]])
			} else {
				x <- c(x[[1:(i[k]-1)]], value[[k]], x[[(i[k]+1):nlyr(x)]])
			}
		}
		x
	}
)



setReplaceMethod("[", c("SpatRaster", "missing", "missing"),
	function(x, i, j, value) {

		nl <- nlyr(x)

		if (is.matrix(value)) {
			if (all(dim(value) == c(ncell(x), nl))) {
				x <- try( setValues(x, value, TRUE, TRUE) )
			} else {
				error(" [,SpatRaster","dimensions of the matrix do not match the SpatRaster")
			}
		} else {
			v <- try( matrix(nrow=ncell(x), ncol=nl) )
			if (! inherits(x, "try-error")) {
				v[] <- value
				x <- try( setValues(x, v, TRUE, TRUE) )
			}
		}
		if (inherits(x, "try-error")) {
			error(" [,SpatRaster", "cannot set values")
		}
		return(x)
	}
)



setReplaceMethod("[", c("SpatRaster","numeric", "missing"),
	function(x, i, j, value) {
		theCall <- sys.call(-1)
		narg <- length(theCall)-length(match.call(call=sys.call(-1)))
		if (narg > 0) { # row
			i <- cellFromRowColCombine(x, i, 1:ncol(x))
		}

		if (!is.null(dim(value))) {
			#x@ptr <- x@ptr$replaceValues(i, value, ncol(value))
			stopifnot(ncol(value) == nlyr(x))
		} else if (length(i) != length(value)) {
			# recycling with warning
			v <- value
			value <- i
			value[] <- v
		}
		if (hasValues(x)) {
			v <- values(x)
		} else {
			v <- matrix(NA, nrow=ncell(x), ncol=nlyr(x))
		}
		v[i,] <- value
		setValues(x, v, TRUE, TRUE)
	}
)




setReplaceMethod("[", c("SpatRaster", "numeric", "numeric"),
	function(x, i, j, value) {
		i <- cellFromRowColCombine(x, i, j)
		x[i] <- value
		x
	}
)



setReplaceMethod("[", c("SpatRaster","missing", "numeric"),
	function(x, i, j, value) {
		i <- cellFromRowColCombine(x, 1:nrow(x), j)
		x[i] <- value
		x
	}
)



setReplaceMethod("[", c("SpatRaster", "logical", "missing"),
	function(x, i, j, value) {
		i <- which(rep_len(i, ncell(x)))
		x[i] <- value
		x
	}
)


setReplaceMethod("[", c("SpatRaster", "SpatRaster", "ANY"),
	function(x, i, j, value) {
		theCall <- sys.call(-1)
		narg <- length(theCall)-length(match.call(call=sys.call(-1)))
		if (narg > 0) { # row
			error(" [,SpatRaster,SpatRaster", "you cannot use a SpatRaster as a row index")
		}
		if (inherits(value, "SpatRaster")) {
			x <- mask(x, i, maskvalues=TRUE)
			cover(x, value)
		} else {
			if (NCOL(value) > 1) {
				error(" [", "cannot use a data.frame with multiple columns")
			}
			value <- unlist(value)
			if (length(value) > 1) {
				v <- values(x)
				i <- as.logical(values(i))
				if (length(value) == sum(i)) {
					v[i] <- value
				} else {
					v[i] <- value[i]
				}	
				values(x) <- v
				x		
			} else {
				mask(x, i, maskvalues=TRUE, updatevalue=value[1])
			}
		}
	}
)


setReplaceMethod("[", c("SpatRaster", "SpatVector", "missing"),
	function(x, i, j, value) {
		theCall <- sys.call(-1)
		narg <- length(theCall)-length(match.call(call=sys.call(-1)))
		if (narg > 0) { # row
			error(" [,SpatRaster,SpatVector", "you cannot use a SpatVector as a row index")
		}
		if (length(value) > 1) {
			value <- rep_len(value, length.out=length(x))
		} 
		rasterize(i, x, field=value, update=TRUE)
	}
)
