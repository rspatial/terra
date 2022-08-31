
setMethod("lines", signature(x="SpatRaster"),
function(x, mx=10000, ...) {
		if(prod(dim(x)) > mx) {
			error("lines", "too many lines (you can increase the value of mx or use as.polygons)")
		}
		v <- as.polygons(x, dissolve=FALSE, values=FALSE)
		lines(v, ...)
	}
)

setMethod("lines", signature(x="SpatVector"),
	function(x, y=NULL, col, lwd=1, lty=1, arrows=FALSE, alpha=1, ...)  {
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
				lines(cbind(a, b), col=col, lwd=lwd, lty=lty, alpha=alpha, ...)
			}
		} else if (grepl("points", gtype)) {
			points(x, col=col, type="l", lwd=lwd, lty=lty, alpha=alpha, ...)
		} else {
			col <- .getCols(n, col, alpha)
			lwd <- rep_len(lwd, n)
			lty <- rep_len(lty, n)
			g <- lapply(x@ptr$get_linesList(), function(i) { names(i)=c("x", "y"); i } )
			for (i in 1:n) {
				plot.xy(g[[i]], type="l", lty=lty[i], col=col[i], lwd=lwd[i], ...)
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
	function(x, col, cex=1, pch=20, alpha=1, ...)  {
		if (nrow(x) == 0) return(invisible(NULL))
		if (missing(col)) col <- "black"
		n <- length(x)
		col <- .getCols(n, col, alpha)
		cex <- rep_len(cex, n)
		pch <- rep_len(pch, n)
		g <- geom(x, df=TRUE)
		if (any(table(g$id) > 1)) {
			g <- split(g, g[,1])
			for (i in 1:n) {
				graphics::points(g[[i]][,3:4], col=col[i], pch=pch[i], cex=cex[i], ...)
			}
		} else {
			graphics::points(g[,3:4], col=col,  pch=pch, cex=cex,...)
		}
	}
)


setMethod("polys", signature(x="SpatVector"),
	function(x, col, border="black", lwd=1, lty=1, alpha=1, ...)  {
		if (nrow(x) == 0) return(invisible(NULL))
		gtype <- geomtype(x)
		if (gtype != "polygons") {
			error("polys", "expecting polygons")
		}
		if (missing(col)) {
			col <- NULL
		}
		cols <- .getCols(length(x), col, alpha)
		out <- list(main_cols=cols)
		out$leg$border <- border
		.plotPolygons(x, out, lwd=lwd, lty=lty, ...)
	}
)

