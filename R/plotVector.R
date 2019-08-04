# Author: Robert J. Hijmans
# Date :  June 2019
# Version 1.0
# License GPL v3

.getCols <- function(n, cols) {
	if (!is.null(cols)) {
		ncols <- length(cols)
		if (ncols > n) {
			steps <- ncols/n
			i <- round(seq(1, ncols, steps))
			cols <- cols[i]				
		} else if (ncols < n) {
			cols <- rep_len(cols, n)
		}
	} 
	cols
}

setMethod("plot", signature(x="SpatVector", y="missing"), 
	function(x, y, col=NULL, xlab="", ylab="", axes=TRUE, add=FALSE, ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		if (couldBeLonLat(x, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} else {
			asp <- 1
		}
		col <- .getCols(size(x), col)
		
		if (gtype == "points") {
			if (is.null(col)) col = "black"
			if (add) {
				points(g[,3], g[,4], col=col, ...)			
			} else {
				plot(g[,3], g[,4], col=col, axes=axes, xlab=xlab, ylab=ylab, asp=asp, ...)
			}
		} else {
			e <- matrix(as.vector(ext(x)), 2)
			if (!add) {
				plot(e, type="n", axes=axes, xlab=xlab, ylab=ylab, asp=asp, ...)
			}
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
					graphics::polypath(a[,3:4], col=col[i], rule = "evenodd", ...)
				}

			} else {
				g <- lapply(g, function(x) split(x, x[,2]))
				#p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4], ...))	
				for (i in 1:length(g)) {
					x <- g[[i]]
					for (j in 1:length(x)) {
						lines(x[[j]][,3:4], col=col[i], ...)
					}
				}
			}
		}
	}
)


setMethod("lines", signature(x="SpatVector"), 
	function(x, col="black", ...)  {
		g <- geom(x)
		gtype <- geomtype(x)
		col <- .getCols(size(x), col)
		if (gtype == "points") {
			graphics::points(g[,3:4], col=col, ...)
		} else {
			g <- split(g, g[,1])
			if (gtype == "polygons") {
				g <- lapply(g, function(x) split(x, x[,c(2,5)]))
			} else {
				g <- lapply(g, function(x) split(x, x[,2]))
			}
			#p <- sapply(g, function(x) lapply(x, function(y) graphics::lines(y[,3:4], ...)))
			for (i in 1:length(g)) {
				x <- g[[i]]
				for (j in 1:length(x)) {
					lines(x[[j]][,3:4], col=col[i], ...)
				}
			}
			
			
		}
	}
)



setMethod("points", signature(x="SpatVector"), 
	function(x, col="black", ...)  {
		col <- .getCols(size(x), col)
		graphics::points(geom(x)[,3:4], col=col, ...)
	}
)

