

.cont.legend <- function(x, ...) {

	cols <- rev(x$leg$cols)
	nc <- length(cols)
	e <- x$leg$ext
	zlim <- x$leg$minmax
	zz <- pretty(zlim, n =(x$leg$levels+1))	
	zz <- zz[zz >= zlim[1] & zz <= zlim[2]]

	if (x$leg$vertical) {
		Y <- seq(e$ymin, e$ymax, length.out=nc+1)
		graphics::rect(e$xmin, Y[-(nc + 1)], e$xmax, Y[-1], col=rev(cols), border=NA, xpd=TRUE)
		ypos <- e$ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dy
		graphics::segments(e$xmin, ypos, e$xmax+e$dx*0.25, ypos, xpd=TRUE)
		text(e$xmax, ypos, formatC(zz, digits=x$leg$digits, format = "f"), pos=4, xpd=TRUE, ...)
	} else {
		X <- seq(e$xmin, e$xmax, length.out=nc+1)
		graphics::rect(X[-(nc + 1)], e$ymin, X[-1], e$ymax, col=rev(cols), border=NA, xpd=TRUE)
		xpos <- e$xmin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dx
		graphics::segments(xpos, e$ymin, xpos, e$ymax+e$dy*0.25, xpd=TRUE)
		text(xpos, e$ymax, formatC(zz, digits=x$leg$digits, format = "f"), pos=1, xpd=TRUE, ...)
	}	
	graphics::rect(e$xmin, e$ymin, e$xmax, e$ymax, border ="black", xpd=TRUE)
}	



.get.leg.coords <- function(x) {

	if (is.null(x$leg$ext)) {
		ext <- x$ext
		xmin <- ext[1]
		xmax <- ext[2]
		ymin <- ext[3]
		ymax <- ext[4]
	} else {
		p <- x$leg$ext
		xmin <- p[1]
		xmax <- p[2]
		ymin <- p[3]
		ymax <- p[4]
		#ymin <- max(ymin, ext["ymin"])
		#ymax <- min(ymax, ext["ymax"])
	}

	leg.shrink <- rep_len(x$leg$shrink,2)
	if (!is.null(x$leg$main)) {
		n <- length(x$leg$main)		
		leg.shrink[2] <- max(x$leg$shrink[2], (.05*n)) 
	}

	yd <- ymax - ymin
	ymin <- ymin + yd * leg.shrink[1]
	ymax <- ymax - yd * leg.shrink[2]
    dx <- xmax - xmin
	dy <- ymax - ymin

	x$leg$ext <- data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, dx=dx, dy=dy)
	x
}

.line.usr <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

.get.leg.extent <- function(x) {
	usr <- graphics::par()$usr
	dx <- graphics::par()$cxy[1] * graphics::par("cex")	
	dy <- graphics::par()$cxy[2] * graphics::par("cex")	
	loc <- x$leg$loc
	p <- NULL
	if (is.character(loc)) {
		if (loc == "right") {
			p <- c(usr[2]+dx, usr[2]+2*dx, usr[3], usr[4])
			x$leg$vertical <- TRUE
		} else if (loc == "left") {
			s <- .line.usr(trunc(par("mar")[2]), 2)
			p <- c(s, s+dx, usr[3], usr[4])
			x$leg$vertical <- TRUE
		} else if (loc == "bottom") {
			s <- .line.usr(trunc(par("mar")[2]), 1)
			p <- c(usr[1], usr[2], s, s+dy)
			x$leg$vertical <- FALSE
		} else {
			stop(loc)
		}
	}
	x$leg$ext <- p
	.get.leg.coords(x)
}

.plot.legend <- function(x) {
	x <- .get.leg.extent(x)
	if (x$leg$type == "continuous") {
		.cont.legend(x)
	}
}


.plotit <- function(x, ...) {

	graphics::par(mar=x$mar)	
	
	plot(x$ext[1:2], x$ext[3:4], type = "n", xlab = "", ylab = "", asp=x$asp, ...)
	rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	
	if (x$leg$legend) {	
		.plot.legend(x)
	}
}	


.as.raster.continuous <- function(x, cols, minmax=NULL, digits, ...) {
		
	out <- list()

	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA

	z <- stats::na.omit(as.vector(Z))
	if (length(z) == 0) stop("no values")
	if (is.null(minmax)) {
		minmax <- range(z)
	}
	interval <- (minmax[2]-minmax[1])/(length(cols)-1)
	breaks <- minmax[1] + interval * (0:(length(cols)-1))
		
	Z[] <- cols[as.integer(cut(Z, breaks, include.lowest=TRUE, right=FALSE))]
	out$raster <- as.raster(Z)

	out$leg <- list()
	out$leg$minmax <- minmax
	out$leg$cols <- cols
	if (missing(digits)) {
		dif <- diff(minmax)
		if (dif == 0) {
			digits = 0;
		} else {
			digits <- max(0, -floor(log10(dif/10)))
		}
	}
	out$leg$digits <- digits
	out$leg$type <- "continuous"
	out$leg$levels <- 5 
	out$leg$shrink <- c(0,0)
	out$leg$main <- NULL
	out$leg$main.cex <-  1	
	
	out
}


getPlotStuff <- function(x, type, cols, maxcell, mar, leg,  ...) {
	out <- list()
	out$mar <- mar
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
	} else {
		out$asp <- 1
	}
	out$ext <- as.vector(ext(x))

	if (type=="classes") {
		ras <- .as.raster.classes(x, cols, ...)
	} else if (type=="range") {
		ras <- .as.raster.range(x, cols, ...)
	} else {
		ras <- .as.raster.continuous(x, cols, ...)
	}
	out$r <- ras$r
	out$leg <- ras$leg
	if (is.na(leg)) {
		out$leg$legend <- FALSE
	} else {
		out$leg$legend <- TRUE
		out$leg$loc <- leg	
	}
	out
}



#setMethod("plot", signature(x="SpatRaster", y="numeric"), 

.plt <- function(x, y=1, maxcell=50000, col=rev(terrain.colors(255)), type="continuous", mar=c(5.1, 4.1, 4.1, 7.1), legend="right", interpolate=FALSE, ...) {

		x <- x[[y]]
		if (!hasValues(x)) {
			stop("SpatRaster has no cell values")
		}
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		x <- getPlotStuff(object, type=type, cols=col, mar=mar, leg=legend, ...)
		x$interpolate <- interpolate

		y <- getPlotStuff(object, type="continuous", cols=rainbow(25), mar=rep(3,4), leg="bottom")
		y$interpolate <- F

		.plotit(x, ...)
	}
#}

