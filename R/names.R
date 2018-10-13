# Author: Robert J. Hijmans
# Date:  October 2008
# Version 0.9
# Licence GPL v3


.uniqueNames <- function(x, sep='.') {
	y <- as.matrix(table(x))
	y <- y[y[,1] > 1, ,drop=F]
	if (nrow(y) > 0) {
		y <- rownames(y)
		for (i in 1:length(y)) {
			j <- which(x==y[i])
			x[j] <- paste(x[j], sep, 1:length(j), sep='')
		}
	}
	x
}



.uniqueNames <- function(x, sep='.') {
	y <- as.matrix(table(x))
	y <- y[y[,1] > 1, ,drop=F]
	if (nrow(y) > 0) {
		y <- rownames(y)
		for (i in 1:length(y)) {
			j <- which(x==y[i])
			x[j] <- paste(x[j], sep, 1:length(j), sep='')
		}
	}
	x
}


.validNames <- function(x, prefix='lyr') {
	x <- trimws(as.character(x))
	x[is.na(x)] <- ""
#	if (.standardnames()) {
		x[x==''] <- prefix
		x <- make.names(x, unique=FALSE)
#	}
	.uniqueNames(x)
}



	
setMethod('names', signature(x='SpatRaster'), 
	function(x) { 
		x@ptr$names
	}
)


setMethod('names<-', signature(x='SpatRaster'), 
	function(x, value)  {
		nl <- nlayer(x)
		if (length(value) != nl) {
			stop('incorrect number of layer names')
		}
		v <- .validNames(value)
		if (!all(v == value)) {
			warning('one or more names were changed to make them valid')
		}
		x@ptr$names <- v
		return(x)
	}
)

