# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3


setMethod("interpolate", signature(object="SpatRaster"), 
	function(object, model, fun=predict, ..., xyNames=c("x", "y"), factors=NULL, const=NULL, na.rm=FALSE, filename="", overwrite=FALSE, wopt=list()) {
	
		out <- rast(object)
		hv <- hasValues(object)
		nms <- c(xyNames, names(object))
		if (length(unique(nms)) != length(nms)) {
			tab <- table(nms)
			stop("duplicate names: ", tab[tab>1])
		}
		nc <- ncol(out)
		testrow <- round(0.5*nrow(object))
		xy <- xyFromCell(out, cellFromRowCol(out, testrow, 1):cellFromRowCol(out, testrow, min(nc, 500)))
		names(xy) <- xyNames
		if (hv) { 
			d <- readValues(object, testrow, 1, 1, nc, TRUE, TRUE)
			xy <- cbind(xy, d)
		}
		r <- .runModel(model, fun, xy, 1, const, (na.rm & hv), ...)
		nl <- ncol(r)
		out <- rast(object, nlyr=nl)
		cn <- colnames(r)
		if (length(cn) == nl) names(out) <- make.names(cn, TRUE)
		
		if (hv) readStart(object)
		b <- writeStart(out, filename, overwrite, wopt)
		for (i in 1:b$n) {
			xy <- xyFromCell(out, cellFromRowCol(out, b$row[i], 1):cellFromRowCol(out, b$row[i]+b$nrows[i]-1, nc))
			names(xy) <- xyNames
			if (hv) { 
				d <- readValues(object, b$row[i], b$nrows[i], 1, nc, TRUE, TRUE)
				xy <- cbind(xy, d)
			}
			v <- .runModel(model, fun, xy, nl, const, (na.rm & hv), ...)
			writeValues(out, v, c(b$row[i], b$nrows[i]))			
		}
		if (hv) readStop(object)
		out <- writeStop(out)
		return(out)
	}
)


