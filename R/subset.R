# Authors: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# License GPL v3



setMethod("subset", signature(x="SpatRaster"), 
function(x, subset, filename="", overwrite=FALSE, wopt=list(), ...) {
	if (is.character(subset)) {
		i <- stats::na.omit(match(subset, names(x)))
		if (length(i)==0) {
			stop("invalid layer names")
		} else if (length(i) < length(subset)) {
			warning("invalid layer names omitted")
		}
		subset <- i
	}

	subset <- as.integer(subset - 1)
	
	opt <- .runOptions(filename, overwrite, wopt)
	x@ptr <- x@ptr$subset(subset, opt)
	show_messages(x, "subset")
	return(x)	
} )


setMethod("$", "SpatRaster",  
	function(x, name) { subset(x, name) } )


setMethod("[[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})


setMethod("[[", c("SpatRaster", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	subset(x, i, ...)
})




setMethod("subset", signature(x="SpatVector"), 
	function(x, subset, drop=FALSE, ...) {
		x <- x[which(as.vector(subset)), , drop=drop]
		show_messages(x, "subset")
	}
)


.subset_cols <- function(x, subset, drop=FALSE) {
	if (is.character(subset)) {
		i <- stats::na.omit(match(subset, names(x)))
	} else {
		i <- subset[subset > 0 & subset <= ncol(x)]
	}
	if (length(i)==0) {
		i = 0
	} 
	if (length(i) < length(subset)) {
		warning("invalid columns omitted")
	}
	x@ptr <- x@ptr$subset_cols(i-1)
	x <- show_messages(x, "subset")
	if (drop) {	# drop geometry
		d <- x@ptr$getDF()
		as.data.frame(d)
	} else {
		x
	}
}

setMethod("$", "SpatVector",  function(x, name) { 
	s <- .subset_cols(x, name, drop=TRUE) 
	s[,1,drop=TRUE]
})


setMethod("[[", c("SpatVector", "numeric", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	.subset_cols(x, i, ..., drop=drop)
})


setMethod("[[", c("SpatVector", "character", "missing"),
function(x, i, j, ... ,drop=FALSE) {
	.subset_cols(x, i, ..., drop=drop)
})


setMethod("[", c("SpatVector", "numeric", "missing"),
function(x, i, j, ... , drop=FALSE) {
	x@ptr <- x@ptr$subset_rows(i-1)
	x <- show_messages(x)
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})


setMethod("[", c("SpatVector", "numeric", "numeric"),
function(x, i, j, ... , drop=FALSE) {
	p <- x@ptr$subset_rows(i-1)
	x@ptr <- p$subset_cols(j-1)	
	x <- show_messages(x)
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})


setMethod("[", c("SpatVector", "missing", "numeric"),
function(x, i, j, ... , drop=FALSE) {
	x@ptr <- x@ptr$subset_cols(j-1)	
	x <- show_messages(x)
	if (drop) {
		as.data.frame(x)
	} else {
		x
	}
})

setMethod("[", c("SpatVector", "missing", "character"),
function(x, i, j, ... , drop=FALSE) {
	m <- match(j, names(x))
	m <- na.omit(m)
	if (length(m) == 0) { 
		m <- 0
	}
	x[,m,drop=drop]
})

setMethod("[", c("SpatVector", "numeric", "character"),
function(x, i, j, ... , drop=FALSE) {
	m <- match(j, names(x))
	m <- na.omit(m)
	if (length(m) == 0) m[,m,drop=drop]
})

