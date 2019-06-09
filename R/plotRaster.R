

.contLegend <- function(p, cols, zlim,  
    digits=0, leg.levels=5, leg.shrink=c(0.1, 0.1), 
	leg.main = NULL, leg.main.cex = 1, ...) {

    nc <- length(cols)
	xmin = p[1]
	xmax = p[2]
	shrk <- (p[4]-p[3])
	leg.shrink <- rep_len(leg.shrink,2)
	if (!is.null(leg.main)) {
		n <- length(leg.main)		
		leg.shrink[2] <- max(leg.shrink[2], (.05*n)) 
	}
	ymin = p[3] + shrk * leg.shrink[1]
	ymax = p[4] - shrk * leg.shrink[2]

    Y <- seq(ymin, ymax, length.out = nc + 1)
    graphics::rect(xmin, Y[-(nc + 1)], xmax, Y[-1], col = rev(cols), border = NA, xpd=TRUE)
    graphics::rect(xmin, ymin, xmax, ymax, border ="black", xpd=TRUE)
	
    dx <- xmax - xmin
	dy <- ymax - ymin
	zz <- pretty(zlim, n = (leg.levels + 1))	
	zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
    ypos <- ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * dy
    graphics::segments(xmin, ypos, xmax + dx * 0.25, ypos, xpd=TRUE)
    text(xmax, ypos, formatC(zz, digits = digits, format = "f"), pos=4, xpd=TRUE, ...)
    if (!is.null(leg.main)) {
		n <- length(leg.main)
		ymax = ymax + 0.05 * dy
        for (i in 1:n) {
			text(x = xmax,
				y = ymax + (n-i) * 0.05 * dy,
				labels = leg.main[i], 
				#adj = c(0.5, 0.5), 
				cex = leg.main.cex,
				xpd=TRUE)
		}
    } 	 	
}


setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function (x, y, cols, maxpixels = 100000, leg.mar=3, 
		leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1,
		useRaster = TRUE, zlim, xlab="", ylab="", ...) {
		
		asp <- 1
		if (couldBeLonLat(x, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} 
		
		if (missing(cols)) cols = rev(grDevices::terrain.colors(25))
		object <- sampleRegular(x[[y]], maxpixels)
		Y <- yFromRow(object, nrow(object):1)
		Z <- t(as.matrix(object, TRUE)[nrow(object):1, , drop = FALSE])
		X <- xFromCol(object, 1:ncol(object))

		if (missing(zlim)) {
			zlim <- range(Z, na.rm=TRUE)
		} else {
			Z[Z < zlim[1]] <- zlim[1]
			Z[Z > zlim[2]] <- zlim[2]
		}
		graphics::plot.new()
		mars <- graphics::par("mar")
		graphics::par(mar = mars + c(0, 0, 0, leg.mar))
		image(X, Y, Z, col = cols, useRaster = useRaster, asp=asp, xlab=xlab, ylab=ylab, ...)
		usr <- graphics::par()$usr
		dx <- graphics::par()$cxy[1] * graphics::par("cex")
		p <- c(usr[2] + 1 * dx, usr[2] + 2 * dx, usr[3], usr[4])
		
		.contLegend(p, cols, zlim, leg.levels=leg.levels, leg.shrink=leg.shrink,
		leg.main=leg.main, leg.main.cex=leg.main.cex)
		
		setHook("before.plot.new", 
			function(...) {
				m <- graphics::par()$mar
				m[4] <- mars[4]
				graphics::par(mar=m)
				setHook("before.plot.new", NULL, action="replace")
			}, action="replace")		
	}
)



setMethod("plot", signature(x="SpatRaster", y="missing"), 
	function(x, y, maxpixels=100000, nc, nr, main, maxnl=16, ...)  {
		nl <- min(nlyr(x), maxnl)
		if (nl == 0) {
			stop('SpatRaster has no cell values')
		}
		if (nl==1) {
			plot(x, 1, ...)
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

		maxpixels=maxpixels/(nl/2)
			
		if (missing(main)) {
			main <- names(x)
		} else {
			main <- rep_len(main, nl)
		}		
			
		for (i in 1:nl) {
			image(x[[i]], main=main[i], ...)
		}
	}
)



setMethod("lines", signature(x="SpatRaster"),
function(x, mx=50000, ...) {
	if(prod(dim(x)) > mx) {
		stop('too many lines')
	}
	v <- as.polygons(x)
	lines(v, ...)
}
)
