

if (!isGeneric("SpatPolygon") ) { setGeneric("SpatPolygon", function(x, ...) standardGeneric("SpatPolygon")) }

setMethod("SpatPolygon", signature(x='missing'), 
	function(...) {
		p <- methods::new('SpatVector')
		p@ptr <- SpatRaster$new()
		return(p)
	}
)


setMethod("SpatPolygon", signature(x='matrix'), 
	function(x, attr=NULL, crs=NA, ...) {
		SpatPolygon( data.frame(x), attr=attr, crs=crs, ...)
	}
)
	
	
setMethod("SpatPolygon", signature(x='data.frame'), 
	function(x, attr=NULL, crs=NA, ...) {

		if (is.null(x$part)) x$part <- 1L
		if (is.null(x$object)) x$object <- 1

		ppp <- SpatPolygons$new()
		x <- split(x, x$object)
		for (i in 1:length(x)) {
			y <- x[[i]]
			pp <- SpatPoly$new()
			if ( any(y$hole > 0) ) {
				ym <- y[y$hole < 1, ]
				z <- split(ym, ym$part)
				for (j in 1:length(z)) {
					p <- SpatPolyPart$new()
					p$set(z[[j]]$x, z[[j]]$y)
					z[[j]] <- p
				}
				yh <- y[y$hole > 0, ]
				zz <- split(yh, yh$part)
				for (j in 1:length(z)) {
					id <- zz[[j]]$hole[1]
					z[[id]]$setHole(zz[[j]]$x, zz[[j]]$y)
				}
				for (j in 1:length(z)) {
					pp$addPart(z[[j]])			
				}
				
			} else {
				z <- split(y, y$part)
				for (j in 1:length(z)) {
					p <- SpatPolyPart$new()
					p$set(z[[j]]$x, z[[j]]$y)
					pp$addPart(p)
				}
			}
			ppp$addPoly(pp)
		}
		
		if (!is.na(crs)) {
			ppp$crs <- crs
		}

		x <- methods::new("SpatVector")
		x@ptr <- ppp
		x
	}
)	


	
	