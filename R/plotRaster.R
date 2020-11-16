

setMethod("plot", signature(x="SpatRaster", y="missing"), 
	function(x, y, maxcell=50000, nc, nr, main, maxnl=16, fun=NULL, ...)  {

		nl <- max(1, min(nlyr(x), maxnl))

		if (nl==1) {
			if (missing(main)) {
				plot(x, 1, maxcell=maxcell, ...)
			} else {
				plot(x, 1, maxcell=maxcell, main=main, ...)
			}
			return(invisible(NULL))
		}
		if (missing(nc)) {
			nc <- ceiling(sqrt(nl))
		} else {
			nc <- max(1, min(nl, round(nc)))
		}
		if (missing(nr)) {
			nr <- ceiling(nl / nc)
		} else {
			nr <- max(1, min(nl, round(nr)))
			nc <- ceiling(nl / nr)
		}
		
		old.par <- graphics::par(no.readonly = TRUE) 
		on.exit(graphics::par(old.par))
		graphics::par(mfrow=c(nr, nc), mar=c(2, 2, 2, 4))

		maxcell=maxcell/(nl/2)
			
		if (missing("main")) {
			main <- names(x)
		} else {
			main <- rep_len(main, nl)	
		}
		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		for (i in 1:nl) {
			#	image(x[[i]], main=main[i], ...)
			plot(x, i, main=main[i], ...)
			if (!is.null(fun)) {
				fun()
			}
		}
	}
)



setMethod("lines", signature(x="SpatRaster"),
function(x, mx=50000, ...) {
	if(prod(dim(x)) > mx) {
		stop("too many lines")
	}
	v <- as.polygons(x)
	lines(v, ...)
}
)
