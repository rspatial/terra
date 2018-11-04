
if (!isGeneric("geomtype")) {setGeneric("geomtype", function(x,...) standardGeneric("geomtype"))}	
if (!isGeneric("geometry")) {setGeneric("geometry", function(x,...) standardGeneric("geometry"))}	


setMethod ('length' , 'SpatLayer', 
	function(x) {
		x@ptr$size()
	}
)

setMethod('geomtype', signature(x='SpatLayer'), 
	function(x, ...){ 
		x@ptr$type()
	}
)	

setMethod('geometry', signature(x='SpatLayer'), 
	function(x, ...){ 
		x@ptr$getGeometry()
	}
)	

setMethod('as.data.frame', signature(x='SpatLayer'), 
	function(x, ...) {
		d <- data.frame(x@ptr$getAttributes(), stringsAsFactors=FALSE)
		names(d) <- x@ptr$names()
		d
	}
)


if (!isGeneric("SpatPolygon") ) { setGeneric("SpatPolygon", function(x, ...) standardGeneric("SpatPolygon")) }

setMethod("SpatPolygon", signature(x='missing'), 
	function(...) {
		p <- methods::new('SpatLayer')
		p@ptr <- SpatLayer$new()
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
			pp <- SpatGeomRings$new()
			if ( any(y$hole > 0) ) {
				ym <- y[y$hole < 1, ]
				z <- split(ym, ym$part)
				for (j in 1:length(z)) {
					p <- SpatGeomRing$new()
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
					p <- SpatGeomRing$new()
					p$set(z[[j]]$x, z[[j]]$y)
					pp$addPart(p)
				}
			}
			ppp$addPoly(pp)
		}
		
		if (!is.na(crs)) {
			ppp$crs <- crs
		}

		x <- methods::new("SpatLayer")
		x@ptr <- ppp
		x
	}
)	


	
	