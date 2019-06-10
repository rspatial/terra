

.plotLegMain <- function(leg.main, xmax, ymax, dy, leg.main.cex) {
    if (!is.null(leg.main)) {
		n <- length(leg.main)
		ymax <- ymax + 0.05 * dy
		for (i in 1:n) {
			text(x=xmax, y=ymax+(n-i)*0.05*dy,
				labels = leg.main[i], cex = leg.main.cex, xpd=TRUE)
		}
	}
}

.getLegCoords <- function(p, ext, leg.shrink, leg.main) {
	xmin <- p[1]
	xmax <- p[2]
	shrk <- (p[4]-p[3])
	leg.shrink <- rep_len(leg.shrink,2)
	if (!is.null(leg.main)) {
		n <- length(leg.main)		
		leg.shrink[2] <- max(leg.shrink[2], (.05*n)) 
	}
	ymin = p[3] + shrk * leg.shrink[1]
	ymax = p[4] - shrk * leg.shrink[2]
	ymin <- max(ymin, ext["ymin"])
	ymax <- min(ymax, ext["ymax"])
    dx <- xmax - xmin
	dy <- ymax - ymin

	data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, dx=dx, dy=dy)
}


.contLegend <- function(p, ext, cols, zlim, digits, leg.levels, leg.shrink, leg.main, leg.main.cex, ...) {

	e <- .getLegCoords(p, ext, leg.shrink, leg.main)
    nc <- length(cols)

    Y <- seq(e$ymin, e$ymax, length.out=nc+1)
    graphics::rect(e$xmin, Y[-(nc + 1)], e$xmax, Y[-1], col=rev(cols), border=NA, xpd=TRUE)
    graphics::rect(e$xmin, e$ymin, e$xmax, e$ymax, border ="black", xpd=TRUE)
	
	zz <- pretty(zlim, n =(leg.levels+1))	
	zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
    ypos <- e$ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dy
    graphics::segments(e$xmin, ypos, e$xmax+e$dx*0.25, ypos, xpd=TRUE)
    text(e$xmax, ypos, formatC(zz, digits=digits, format = "f"), pos=4, xpd=TRUE, ...)
	.plotLegMain(leg.main, e$xmax, e$ymax, e$dy, leg.main.cex)
}

.sampleColors <- function(cols, n) {
	if (length(cols) != n) {
		if (n==1) {
			cols <- cols[round(length(cols)/2)]
		} else if (n==2) {
			cols <- c(cols[1], cols[length(cols)])
		} else {
			colstep <- (length(cols)-1) / (n-1)
			i <- round(seq(1, length(cols), colstep))
			cols <- cols[i]
		}
	}
	cols
}



.fewClassLegend <- function(p, ext, u, cols, digits, leg.shrink, leg.main, leg.main.cex, ...) {
	u <- sort(u)
	n <- length(u)
	e <- .getLegCoords(p, ext, leg.shrink, leg.main)
	step <- e$dy/20
    Y <- e$ymax - 0:(n-1) * (step *1.5)
	#to put in the middle
	#mid <- Y[trunc((length(Y)+1)/2)]
	#shift <- max(0, mid - ((e$ymax - e$ymin) / 2))
	#Y <- Y - shift
	
	for (i in 1:n) {
		graphics::rect(e$xmin, Y[i], e$xmax, Y[i]-step, col=cols[i], border="black", xpd=TRUE)
	}
    text(e$xmax, Y-0.5*step, formatC(u, digits=digits, format = "f"), pos=4, xpd=TRUE, ...)
	.plotLegMain(leg.main, e$xmax, e$ymax, e$dy, leg.main.cex)	
}


setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function (x, y, cols, maxpixels = 100000, leg.mar=3, 
		leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1,
		digits, useRaster = TRUE, zlim, xlab="", ylab="", ...) {
		
		asp <- 1
		if (couldBeLonLat(x, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} 
		
		if (missing(cols)) cols = rev(grDevices::terrain.colors(25))
		object <- sampleRegular(x[[y]], maxpixels)
		xex <- as.vector(ext(object))
		
		Y <- yFromRow(object, nrow(object):1)
		Z <- t(as.matrix(object, TRUE)[nrow(object):1, , drop = FALSE])
		X <- xFromCol(object, 1:ncol(object))

		if (missing(zlim)) {
			zlim <- range(Z, na.rm=TRUE)
		} else {
			zlim <- sort(zlim)
			Z[Z < zlim[1]] <- zlim[1]
			Z[Z > zlim[2]] <- zlim[2]
		}
		dv <- dev.list()
		if (!is.null(dv)) {
			# excluding pdf, png, ....
			# what are other others ones to include beside windows?
			if (names(dv[length(dv)]) %in% c("windows")){
				graphics::plot.new()
			}
		}
#		graphics::plot.new()
		mars <- graphics::par("mar")

		if (missing(digits)) {
			dif <- diff(zlim)
			if (dif == 0) {
				digits = 0;
			} else {
				digits <- max(0, -floor(log10(dif/10)))
			}
		}
		u <- unique(na.omit(as.vector(Z)))
		if (length(u) < 10) {
			cols <- .sampleColors(cols, length(u))
		}
		
		graphics::par(mar = mars + c(0, 0, 0, leg.mar))
		image(X, Y, Z, col = cols, useRaster = useRaster, asp=asp, xlab=xlab, ylab=ylab, ...)

		usr <- graphics::par()$usr
		dx <- graphics::par()$cxy[1] * graphics::par("cex")
		p <- c(usr[2] + 1 * dx, usr[2] + 2 * dx, usr[3], usr[4])
	
		u <- unique(as.vector(Z), na.rm=TRUE)
		if (length(u) < 10) {
			.fewClassLegend(p, xex, u, cols, digits, leg.shrink, leg.main, leg.main.cex)		
		} else {
			.contLegend(p, xex, cols, zlim, digits, leg.levels, leg.shrink, leg.main, leg.main.cex)
		}
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
