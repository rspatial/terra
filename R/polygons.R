
	


.spatPolygons <- function(x, ..., attr=NULL, crs=NA) {
	Part <- function(xy) {
		p <- SpatPolyPart$new()
		p$set(xy[,1], xy[,2])
		p
	}

	Parts <- function(x) {
		pp <- SpatPoly$new()
		if (length(x) == 1) {
			pp$addPart(x)	
		} else {
			for (i in 1:length(x)) {
				pp$addPart(x[[i]])
			}
		}
		pp
	}
	
	x <- c(list(x), list(...))
	y <- rapply(x, Part, how='replace')
	z <- lapply(y, Parts)

	ppp <-  SpatPolygons$new()
	lapply(z, function(i) ppp$addPoly(i))

	if (!is.na(crs)) {
		ppp$crs <- crs
	}
	if (!is.null(attr)) {
		ppp$attr <- attr[,1]
	}
	
	x <- methods::new("SpatPolygons")
	x@ptr <- ppp
	x
}


setMethod ('length' , 'SpatPolygons', 
	function(x) {
		x@ptr$size()
	}
)

setMethod('ext', signature(x='SpatPolygons'), 
	function(x, ...){ 
		e <- methods::new('SpatExtent')
		e@ptr <- x@ptr$extent
		return(e)
	}
)	



setMethod ('show' , 'SpatPolygons', 
	function(object) {
		
		cat('class       :' , class(object), '\n')

		d <- length(object)
		cat('geometries  : ', d, ' \n', sep="" ) 
		e <- as.vector(ext(object))
		cat('extent      : ' , e[1], ', ', e[2], ', ', e[3], ', ', e[4], '  (xmin, xmax, ymin, ymax)\n', sep="")
		crs <- crs(object)
		cat('coord. ref. :' , crs(object), '\n')
		
	}
)


