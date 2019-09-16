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

.plotPolygons <- function(x, cols, ...) {
	g <- geom(x)
	g <- split(g, g[,1])
	g <- lapply(g, function(y) split(y, y[,2]))
	for (i in 1:length(g)) {
		gg <- g[[i]]
		for (j in 1:length(gg)) {
			a <- gg[[j]]
			if (any(a[,5] > 0)) {
				a <- split(a, a[,5]) 
				a <- lapply(a, function(i) rbind(i, NA))
				a <- do.call(rbind, a )
				a <- a[-nrow(a), ]
				# g[[i]][[1]] <- a 
			}
			graphics::polypath(a[,3:4], col=cols[i], rule = "evenodd", ...)
		}
	}
}

.plotLines <- function(x, cols, ...) {
	g <- geom(x)
	g <- split(g, g[,1])
	g <- lapply(g, function(x) split(x, x[,2]))
	#p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4], ...))	
	for (i in 1:length(g)) {
		x <- g[[i]]
		for (j in 1:length(x)) {
			lines(x[[j]][,3:4], col=col[i], ...)
		}
	}
}


setMethod("plot", signature(x="SpatVector", y="missing"), 
	function(x, y, col=NULL, xlab="", ylab="", axes=TRUE, add=FALSE, ...)  {
		gtype <- geomtype(x)
		if (couldBeLonLat(x, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} else {
			asp <- 1
		}
		col <- .getCols(size(x), col)
		
		if (gtype == "points") {
			if (is.null(col)) col = "black"
			g <- geom(x)
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
			if (gtype == "polygons") {
				.plotPolygons(x, col, ...)
			} else {
				.plotLines(x, col, ...)
			}
		}
	}
)


.setPAR <- function(leg.mar=NULL, leg.ext=NULL, leg.hor=NULL, ...) {
	if (missing(leg.mar)) {
		if (missing(leg.ext)) {
			leg.mar=4
		} else {
			leg.mar=0
		}
	}
	leg.hor <- FALSE
	if (leg.hor) {
		graphics::par(mar=.getMar(c(leg.mar, 0, 0, 0)))
	} else {
		graphics::par(mar=.getMar(c(0, 0, 0, leg.mar)))
	}		
}

.legCoords <- function(object, leg.ext=NULL, leg.shrink=c(0,0), leg.main=NULL, ...) {
	usr <- graphics::par()$usr
	dx <- graphics::par()$cxy[1] * graphics::par("cex")	

	if (!is.null(leg.ext)) {
		xex <- as.vector(leg.ext)
		leg.ext <- .getLegCoords(NULL, xex, leg.shrink, leg.main)
	} else {
		p <- c(usr[2]+dx, usr[2]+2*dx, usr[3], usr[4])
		xex <- as.vector(ext(object))
		leg.ext <- .getLegCoords(p, xex, leg.shrink, leg.main)
	} 
	leg.ext
}


setMethod("plot", signature(x="SpatVector", y="character"), 
	function(x, y, col=topo.colors(100), xlab="", ylab="", axes=TRUE, add=FALSE, leg.ext = NULL, ...)  {
		v <- unlist(x[, y, drop=TRUE], use.names=FALSE)
		uv <- unique(v)
		if (is.null(col)) {
			col 
		}
		ucols <- .getCols(length(uv), col)
		i <- match(v, uv)
		cols <- ucols[i]
		
		old.par <- graphics::par(no.readonly = TRUE) 
		on.exit(graphics::par(old.par))
		.setPAR(...)
		plot(x, col=cols, ...)
		n <- ifelse(is.null(leg.ext), 20, length(uv))
		leg.ext <- .legCoords(x, ...)
		.factorLegend(leg.ext, 1:length(uv), ucols, uv, n)
	}
)



setMethod("lines", signature(x="SpatVector"), 
	function(x, col=NULL, ...)  {
		if (is.null(col)) col <- "black"
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
	function(x, col=NULL, ...)  {
		col <- .getCols(size(x), col)
		if (is.null(col)) col <- "black"
		graphics::points(geom(x)[,3:4], col=col, ...)
	}
)

