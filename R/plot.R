# Author: Robert J. Hijmans
# Date :  June 2008
# Version 1.0
# License GPL v3

#setMethod("plot", signature(x="SpatRaster", y="numeric"), 
#	function(x, y, maxpixels=100000, xlab="", ylab="", ...)  {
#		y <- as.integer(y[1])
#		stopifnot(y>0 && y<=nlyr(x))
#		x <- x[[y]]
#		stopifnot(.hasValues(x));
#		x <- sampleRegular(x, maxpixels)
#		m <- matrix(values(x), nrow=nrow(x), byrow=TRUE)	
#		require(lattice)
#		lattice::levelplot(t(m[nrow(m):1, , drop=FALSE]), xlab=xlab, ylab=ylab, ...)
#	}
#)

setMethod("plot", signature(x="SpatRaster", y="missing"), 
	function(x, y, maxpixels=100000, xlab="", ylab="", ...)  {
		plot(x, 1, maxpixels=maxpixels, xlab=xlab, ylab=ylab, ...)
	}
)

setMethod("plot", signature(x="SpatVector", y="missing"), 
	function(x, y, xlab="", ylab="", ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		if (gtype == "points") {
			plot(g[,3], g[,4], xlab=xlab, ylab=ylab, ...)
		} else {
			e <- matrix(as.vector(ext(x)), 2)
			plot(e, type="n", xlab=xlab, ylab=ylab, ...)
			g <- split(g, g[,1])
			if (gtype == "polygons") {
				g <- lapply(g, function(x) split(x, x[,2]))
				for (i in 1:length(g)) {
					a <- g[[i]][[1]]
					if (any(a[,5] > 0)) {
						a <- split(a,a[,5]) 
						a <- lapply(a, function(i) rbind(i, NA))
						a <- do.call(rbind, a )
						a <- a[-nrow(a), ]
						g[[i]][[1]] <- a
					}
					polypath(a[,3:4], rule = "evenodd", ...)
				}

			} else {
				g <- lapply(g, function(x) split(x, x[,2]))
				p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4], ...)))			
			}
		}
	}
)


setMethod("lines", signature(x="SpatVector"), 
	function(x, ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		if (gtype == "points") {
			points(g[,3:4], ...)
		} else {
			g <- split(g, g[,1])
			g <- lapply(g, function(x) split(x, x[,2]))
			p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4], ...)))
		}
	}
)

setMethod("points", signature(x="SpatVector"), 
	function(x, ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		points(g[,3:4], ...)
	}
)


setMethod("lines", signature(x="SpatExtent"), 
	function(x, ...)  {
		e <- as.vector(x)
		p <- rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3]), c(e[1],e[3]))
		lines(p, ...)		
	}
)

