

setMethod("$<-", "SpatRaster",  
	function(x, name, value) { 
		if (inherits(value, "SpatRaster")) {
			value <- value[[1]]
			names(value) <- name
		} else if (!is.null(value)) {
			stop("the replacement value should be a SpatRaster")
		}
	
		i <- which(name == names(x))[1]
		if (is.na(i)) {
			return(c(x, value))
		} else {
			if (i == 1) {
				x <- c(value, x[[2:nlyr(x)]])
			} else if (i == nlyr(x)) {
				x <- c(x[[1:(nlyr(x)-1)]], value)			
			} else {
				x <- c(x[[1:(i-1)]], value, x[[(i+1):nlyr(x)]])
			}
			return(x)
		} 
	}
)


setReplaceMethod("[[", c("SpatRaster", "character", "missing"),
	function(x, i, j, value) {
		if (nlyr(value) != length(i)) {
			stop("length of names must be equal to the number of layers")
		}
		names(value) <- i
		for (k in 1:length(i)) {
			eval(parse(text = paste("`$`(x, ", i[k], ") <- value[[k]]")))
		}
		x
	}
)

setReplaceMethod("[[", c("SpatRaster", "numeric", "missing"),
	function(x, i, j, value) {
		if (nlyr(value) != length(i)) {
			stop("length of indices must be equal to the number of layers")
		}
		if (any(i<1) | any(i > nlyr(x))) {
			stop("indices must be between 1 and the number of layers")
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
		x <- rast(x)
		
		if (is.matrix(value)) {
			if (all(dim(value) == c(ncell(x), nl))) {
				e <- try( values(x) <- value)
			} else {
				stop("dimensions of the matrix do not match the SpatRaster")
			}
		} else {
			v <- try( matrix(nrow=ncell(x), ncol=nl) )
			if (! inherits(x, "try-error")) {
				v[] <- value
				e <- try(values(x) <- v)
			}
		}
		if (inherits(e, "try-error")) {
			stop("cannot set values")
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
		v <- values(x)
		v[i,] <- value
		values(x) <- v
		x
	}
)


setReplaceMethod("[", c("SpatRaster", "numeric", "numeric"),
	function(x, i, j, value) {
		i <- cellFromRowColCombine(x, i, j)
		v <- values(x)
		v[i,] <- value
		values(x) <- v
		x
	}
)	


setReplaceMethod("[", c("SpatRaster","missing", "numeric"),
	function(x, i, j, value) {
		i <- cellFromRowColCombine(x, 1:nrow(x), j)
		v <- values(x)
		v[i,] <- value
		values(x) <- v
		x
	}
)



setReplaceMethod("[", c("SpatRaster", "logical", "missing"),
	function(x, i, j, value) {
		v <- values(x)
		v[i,] <- value
		values(x) <- v
		x
	}
)	

