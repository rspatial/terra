# Author: Robert J. Hijmans
# Date:  October 2018
# Version 1.0
# License GPL v3


.rast_replace <- function(x, name, value, caller="$<-") {
	if (inherits(value, "SpatRaster")) {
		value <- value[[1]]
		names(value) <- name
	} else if (!is.null(value)) {
		y <- rast(x, nlyrs=1)
		test <- try(values(y) <- value, silent=TRUE)
		if (inherits(test, "try-error")) {
			error(caller, "the replacement value is not valid")
		}
		value <- y
		names(value) <- name
	}
	i <- which(name == names(x))[1]
	if (is.null(value)) {
		if (is.na(i)) {
			return(x)
		} else {
			return(subset(x, -i, NSE=FALSE))
		}
	}

	if (is.na(i)) {
		if (hasValues(x)) {
			c(x, value)
		} else if (hasValues(value)) {
			value
		} else {
			c(x, value)			
		}
	} else if (nlyr(x) == 1) {
		value$deepcopy()
	} else if (i == 1) {
		c(value, x[[2:nlyr(x)]])
	} else if (i == nlyr(x)) {
		c(x[[1:(nlyr(x)-1)]], value)
	} else {
		c(x[[1:(i-1)]], value, x[[(i+1):nlyr(x)]])
	}
}

setMethod("$<-", "SpatRaster",
	function(x, name, value) {
		.rast_replace(x, name, value, "`$<-`")
	}
)


setReplaceMethod("[[", c("SpatRaster", "character"),
	function(x, i, value) {
		if (inherits(value, "numeric")) {
			r <- rast(x, nlyr=length(i))
			value <- init(r, value)
			names(value) <- i
		} else if (inherits(value, "SpatRaster")) {
			if (nlyr(value) != length(i)) {
				error("`[[`", "length of names must be equal to the number of layers")
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
			x <- .rast_replace(x, i[k], value[[k]], "`[[<-`")
#			eval(parse(text = paste0("x$", i[k], " <- value[[k]]")))
		}
		x
	}
)

setReplaceMethod("[[", c("SpatRaster", "numeric"),
	function(x, i, value) {
		if (inherits(value, "numeric")) {
			r <- rast(x, nlyr=length(i))
			value <- init(r, value)
		} else if (!inherits(value, "SpatRaster")) {
			error("`[[<-`", "Expected a SpatRaster or numeric as replacement value")
		}
		if (nlyr(value) < length(i)) {
			if (nlyr(value) > 1) {
				j <- rep_len(1:nlyr(value), length(i))
				value <- value[[j]]
			}
		} else if (nlyr(value) > length(i)) {
			error("`[[`", "length of indices must be <= the number of layers")
		}
		if (any(i<1) | any(i > nlyr(x))) {
			error("`[[`", "indices must be between 1 and the number of layers")
		}
		if (nlyr(x) == 1) {
			compareGeom(x, value, crs=FALSE, warncrs=TRUE)
			return(value)
		}
		if (nlyr(value) == 1) {
			for (k in 1:length(i)) {
				if (i[k] == 1) {
					x <- c(value, x[[2:nlyr(x)]])
				} else if (i[k] == nlyr(x)) {
					x <- c(x[[1:(nlyr(x)-1)]], value)
				} else {
					x <- c(x[[1:(i[k]-1)]], value, x[[(i[k]+1):nlyr(x)]])
				}
			}
		} else {
			for (k in 1:length(i)) {
				if (i[k] == 1) {
					x <- c(value[[k]], x[[2:nlyr(x)]])
				} else if (i[k] == nlyr(x)) {
					x <- c(x[[1:(nlyr(x)-1)]], value[[k]])
				} else {
					x <- c(x[[1:(i[k]-1)]], value[[k]], x[[(i[k]+1):nlyr(x)]])
				}
			}
		}
		#g <- gc()
		x
	}
)


