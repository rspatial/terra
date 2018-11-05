

if (!isGeneric("vect") ) { setGeneric("vect", function(x, ...) standardGeneric("vect")) }

setMethod("vect", signature(x='missing'), 
	function(...) {
		p <- methods::new('SpatLayer')
		p@ptr <- SpatLayer$new()
		.messages(p)
		return(p)
	}
)

setMethod("vect", signature(x='character'), 
	function(x, ...) {
		p <- methods::new('SpatLayer')
		p@ptr <- SpatLayer$new()
		p@ptr$read(x)
		.messages(p)
		p
	}
)



setMethod("vect", signature(x='matrix'), 
	function(x, type='points', atts=NULL, crs=NA, ...) {
		type <- tolower(type)
		stopifnot(type %in% c("points", "lines", "polygons"))
		
		p <- methods::new('SpatLayer')
		p@ptr <- SpatLayer$new()
		nr <- nrow(x)
		if (nr == 0) {
			return(p)
		}
		
		if (ncol(x) == 2) {
			p@ptr$setGeometry(type, 1:nr, rep(1, nr), x[,1], x[,2], rep(FALSE, nr))
		} else if (ncol(x) == 4) {
			p@ptr$setGeometry(type, x[,1], x[,2], x[,3], x[,4], rep(FALSE, nr))		
		} else if (ncol(x) == 5) {
			p@ptr$setGeometry(type, x[,1], x[,2], x[,3], x[,4], x[,5])
		} else {
			stop("not an appropriate matrix")
		}
		if (!is.null(atts)) {
			p@ptr$setAttributes(p)
		}
		if (!is.na(crs)) {
			p@ptr$setCrs(crs)
		}
		.messages(p)
		p
	}
)


setMethod("vect", signature(x='data.frame'), 
	function(x, type='points', atts=NULL, crs=NA, ...) {
		x <- as.matrix(x)
		vect(x, type=type, atts=atts, crs=crs, ...)
	}
)
