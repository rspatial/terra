
setMethod("lines", signature(x="SpatRaster"),
function(x, mx=10000, ...) {
		reset.clip()
		if(prod(dim(x)) > mx) {
			error("lines", "too many lines (you can increase the value of mx or use as.polygons)")
		}
		v <- as.polygons(x[[1]], dissolve=FALSE, values=FALSE)
		lines(v, ...)
	}
)


setMethod("points", signature(x="SpatRaster"),
function(x, ...) {
		reset.clip()
		p <- as.points(x[[1]])
		points(p, ...)
	}
)

setMethod("polys", signature(x="SpatRaster"),
function(x, mx=10000, dissolve=TRUE, ...) {
		reset.clip()
		if(prod(dim(x)) > mx) {
			error("lines", "too many lines (you can increase the value of mx or use as.polygons)")
		}
		p <- as.polygons(x[[1]], dissolve=dissolve)
		polys(p, ...)
	}
)


setMethod("lines", signature(x="SpatVector"),
	function(x, y=NULL, col, lwd=1, lty=1, arrows=FALSE, alpha=1, ...)  {
		reset.clip()
		n <- nrow(x)
		if (n == 0) return(invisible(NULL))
		gtype <- geomtype(x)
		if (missing(col)) col <- "black"
		if (!is.null(y)) {
			stopifnot(inherits(y, "SpatVector"))
			ytype <- geomtype(y)
			if ((ytype != "points") || (gtype != "points")) {
				error("lines", "when supplying two SpatVectors, both must have point geometry")
			}
			stopifnot(nrow(x) == nrow(y))
			p1 <- geom(x)[, c("x", "y"), drop=FALSE]
			p2 <- geom(y)[, c("x", "y"), drop=FALSE]
			if (arrows) {
				arrows(p1[,1], p1[,2], p2[,1], p2[,2], col=col, lwd=lwd, lty=lty, ...)
			} else {
				a <- as.vector(t(cbind(p1[,1], p2[,1], NA)))
				b <- as.vector(t(cbind(p1[,2], p2[,2], NA)))
				lines(cbind(a, b), col=col, lwd=lwd, lty=lty, ...)
			}
		} else {
			if (gtype != "polygons") {
				x <- as.lines(x) 
			}
			if ((length(col) == 1) && (length(lty)==1) && (length(lwd)==1)) {
				col <- .getCols(1, col, alpha)
				g <- x@cpp$linesNA()
				names(g) <- c("x", "y")
				graphics::plot.xy(g, type="l", lty=lty, col=col, lwd=lwd, ...)
			} else {
				col <- .getCols(n, col, alpha)
				lwd <- rep_len(lwd, n)
				lty <- rep_len(lty, n)
#				g <- lapply(x@cpp$linesList(), function(i) { names(i)=c("x", "y"); i } )
				g <- x@cpp$linesList()
				for (i in 1:n) {
					if (length(g[[i]]) > 0) {
						names(g[[i]]) = c("x", "y")
						graphics::plot.xy(g[[i]], type="l", lty=lty[i], col=col[i], lwd=lwd[i])
					}
				}
			}
			#g <- geom(x, df=TRUE)
			#g <- split(g, g[,1])
			#if (gtype == "polygons") {
			#	g <- lapply(g, function(x) split(x, x[,c(2,5)]))
			#} else {
			#	g <- lapply(g, function(x) split(x, x[,2]))
			#}
			#for (i in 1:length(g)) {
			#	for (j in 1:length(g[[i]])) {
			#		lines(g[[i]][[j]][,3:4], col=col[i], lwd=lwd[i], lty=lty[i], ...)
			#	}
			#}
		}
	}
)


setMethod("points", signature(x="SpatVector"),
	function(x, col, cex=0.7, pch=16, alpha=1, ...)  {
		reset.clip()
		n <- length(x)
		if (n == 0) return(invisible(NULL))
		if (missing(col)) col <- "black"
		if ((length(col) == 1) && (length(cex)==1) && (length(pch)==1)) {
			col <- .getCols(1, col, alpha)
			#graphics::points(g[,3:4], col=col,  pch=pch, cex=cex,...)
			g <- crds(x)
			graphics::plot.xy(list(x=g[,1], y=g[,2]), type="p", pch=pch, col=col, cex=cex, ...)
		} else {
			col <- .getCols(n, col, alpha)
			cex <- rep_len(cex, n)
			pch <- rep_len(pch, n)
			g <- geom(x, df=TRUE)
			if (nrow(g) > g[nrow(g), 1]) {
				g <- split(g[,3:4], g[,1])
				for (i in 1:n) {
					#graphics::points(g[[i]], col=col[i], pch=pch[i], cex=cex[i], ...)
					graphics::plot.xy(list(x=g[[i]][,1], y=g[[i]][,2]), type="p", pch=pch[i], col=col[i], cex=cex[i], ...)
				}
			} else {
				graphics::plot.xy(list(x=g[,3], y=g[,4]), type="p", pch=pch, col=col, cex=cex, ...)
			}
		}
	}
)


setMethod("polys", signature(x="SpatVector"),
	function(x, col, border="black", lwd=1, lty=1, alpha=1, ...)  {
		reset.clip()
		gtype <- geomtype(x)
		if (gtype != "polygons") {
			error("polys", "expecting polygons")
		}
		if (missing(col)) {
			col <- NULL
		} else if (length(col) > 1) {
			col <- .getCols(length(x), col, alpha)
		}
		out <- list(main_cols=col)
		out$leg$border <- border
		p <- .plotPolygons(x, out, lwd=lwd, lty=lty, ...)
	}
)

		
setMethod("points", signature(x="sf"),
function(x, ...) {
		points(vect(x), ...)
	}
)				
		
setMethod("lines", signature(x="sf"),
function(x, ...) {
		lines(vect(x), ...)
	}
)

setMethod("polys", signature(x="sf"),
function(x, ...) {
		polys(vect(x), ...)
	}
)
