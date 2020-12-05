setMethod("lines", signature(x="SpatRaster"),
function(x, mx=50000, ...) {
		if(prod(dim(x)) > mx) {
			stop("too many lines (you can increase the value of mx)")
		}
		v <- as.polygons(x)
		lines(v, ...)
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

setMethod("polys", signature(x="SpatVector"), 
	function(x, col, border="black", ...)  {
		gtype <- geomtype(x)
		if (gtype != "polygons") {
			stop("expecting polygons")
		}
		if (missing(col)) {
			col <- NULL
		}
		cols <- .getCols(size(x), col)
		.plotPolygons(x, cols=cols, border=border, ...)
	}
)
	
