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

.plotPolygons <- function(x, cols, border=NULL, ...) {
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
			graphics::polypath(a[,3:4], col=cols[i], rule = "evenodd", border=border, ...)
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
			lines(x[[j]][,3:4], col=cols[i])
		}
	}
}

.vplot <- function(x, y, col, axes=TRUE, add=FALSE, border="black", xlab="", ylab="", asp=NULL, ...) {
	gtype <- geomtype(x)
	if (is.null(asp)) {
		if (isLonLat(x, perhaps=TRUE, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} else {
			asp <- 1
		}
	}
	if (gtype == "points") {
		if (missing(col)) col = "black"
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
			.plotPolygons(x, col, border=border, ...)
		} else {
			if (is.null(col)) col = rep("black", size(x))
			.plotLines(x, col, ...)
		}
	}
}

setMethod("plot", signature(x="SpatVector", y="missing"), 
	function(x, y, col, axes=TRUE, add=FALSE, ...)  {
		if (missing(col)) col <- ifelse(geomtype(x) == "points", "black", NULL)
		col <- .getCols(size(x), col)
		.vplot(x, y, col=col, axes=axes, add=add, ...)
	}
)


#.setPAR <- function(leg.mar=NULL, leg.ext=NULL, leg.hor=NULL, ...) {
#	if (missing(leg.mar)) {
#		if (missing(leg.ext)) {
#			leg.mar=4
#		} else {
#			leg.mar=0
#		}
#	}
#	leg.hor <- FALSE
#	if (leg.hor) {
#		graphics::par(mar=.getMar(c(leg.mar, 0, 0, 0)))
#	} else {
#		graphics::par(mar=.getMar(c(0, 0, 0, leg.mar)))
#	}		
#}

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


# note: add=T is not working
# signature(x="SpatVector", y="missing"),
# should be calling this one, and not the other way around? 

setMethod("plot", signature(x="SpatVector", y="character"), 
	function(x, y, col, type, mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, add=FALSE, ...)  {
		
		#old.par <- graphics::par(no.readonly = TRUE) 
		#on.exit(graphics::par(old.par))
		#.setPAR(...)
		
		if (missing(col)) {
			col <- topo.colors(100)
		}

		if (!add) graphics::par(mar=mar)

		v <- unlist(x[, y, drop=TRUE], use.names=FALSE)
		uv <- unique(v)
		if (missing(type)) {
			if (!is.numeric(uv) | length(uv) < 10) {
				type <- "classes"
			} else {
				type <- "continuous"
			}
		}
		if (type == "classes") {
			ucols <- .getCols(length(uv), col)
			uv <- sort(uv)
			i <- match(v, uv)
			cols <- ucols[i]
		} else {
			brks <- seq(min(v, na.rm=TRUE), max(v, na.rm=TRUE), length.out = length(col))
			grps <- cut(v, breaks = brks, include.lowest = TRUE)
			cols=col[grps]
		}

		plot(x, col=cols, add=add, ...)
		leg.ext <- NULL
		n <- ifelse(is.null(leg.ext), 20, length(uv))
		leg.ext <- .legCoords(x, ...)
		if (type == "classes") {
			.factorLegend(leg.ext, 1:length(uv), ucols, uv, n)
		} else {
			zlim <- range(uv, na.rm=TRUE)
			if (missing(digits)) {
				dif <- diff(zlim)
				if (dif == 0) {
					digits = 0;
				} else {
					digits <- max(0, -floor(log10(dif/10)))
				}
			}
			.contLegend(leg.ext, col, zlim, digits, leg.levels=5, ...)	
		}
	}
)



setMethod("lines", signature(x="SpatVector"), 
	function(x, col, ...)  {
		if (missing(col)) col <- "black"
		g <- geom(x)
		gtype <- geomtype(x)
		if (missing(col)) col <- NULL
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
	function(x, col, ...)  {
		if (missing(col)) col <- "black"
		col <- .getCols(size(x), col)
		graphics::points(geom(x)[,3:4], col=col, ...)
	}
)

