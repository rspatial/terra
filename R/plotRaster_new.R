

.cont.legend <- function(x, ...) {

	cols <- rev(x$leg$cols)
	nc <- length(cols)
	e <- x$leg$ext
	zlim <- x$leg$minmax
	zz <- x$leg$at
	if (is.null(zz)) {
		zz <- pretty(zlim, n =(x$leg$levels+1))	
		zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
	}
	
	if (x$leg$loc %in% c("left", "right")) {
		Y <- seq(e$ymin, e$ymax, length.out=nc+1)
		graphics::rect(e$xmin, Y[-(nc + 1)], e$xmax, Y[-1], col=rev(cols), border=NA, xpd=TRUE)
		ypos <- e$ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dy
		if (x$leg$loc == "right") {
			graphics::segments(e$xmin, ypos, e$xmax+e$dx*0.25, ypos, xpd=TRUE)
			text(e$xmax, ypos, formatC(zz, digits=x$leg$digits, format = "f"), pos=4, xpd=TRUE, ...)
		} else {
			graphics::segments(e$xmin-e$dx*0.25, ypos, e$xmax, ypos, xpd=TRUE)
			text(e$xmin, ypos, formatC(zz, digits=x$leg$digits, format = "f"), pos=2, xpd=TRUE, ...)
		}
	} else {
		X <- seq(e$xmin, e$xmax, length.out=nc+1)
		graphics::rect(X[-(nc + 1)], e$ymin, X[-1], e$ymax, col=rev(cols), border=NA, xpd=TRUE)
		xpos <- e$xmin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dx
		if (x$leg$loc == "bottom") {
			graphics::segments(xpos, e$ymin-e$dy*0.25, xpos, e$ymax, xpd=TRUE)
			text(xpos, e$ymin, formatC(zz, digits=x$leg$digits, format = "f"), pos=1, xpd=TRUE, ...)
		} else {
			graphics::segments(xpos, e$ymin, xpos, e$ymax+e$dy*0.25, xpd=TRUE)
			text(xpos, e$ymax+e$dy*0.25, formatC(zz, digits=x$leg$digits, format = "f"), pos=3, xpd=TRUE, ...)
		}
	}	
	graphics::rect(e$xmin, e$ymin, e$xmax, e$ymax, border ="black", xpd=TRUE)
}	

#.plt(r, leg="left", mar=c(3,6,3,3))


.get.leg.coords <- function(x) {

	if (is.null(x$leg$ext)) {
		ext <- unlist(x$ext)
		xmin <- x$ext[1]
		xmax <- x$ext[2]
		ymin <- x$ext[3]
		ymax <- x$ext[4]
	} else {
		p <- unlist(x$leg$ext)
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
	lh <- par("cin")[2] * par("cex") * par("lheight")
	x_off <- diff(grconvertX(c(0, lh), "inches", "npc"))
	y_off <- diff(grconvertY(c(0, lh), "inches", "npc"))
	if (side == 1) {
		grconvertY(-line * y_off, "npc", "user")
	} else if (side ==2) {
		grconvertX(-line * x_off, "npc", "user")
	} else if (side ==3) {
		grconvertY(1 + line * y_off, "npc", "user")
	} else {
		grconvertX(1 + line * x_off, "npc", "user")
	}
}

.get.leg.extent <- function(x) {
	usr <- graphics::par()$usr
	dxy <- graphics::par()$cxy * graphics::par("cex")	
	loc <- x$leg$loc
	p <- NULL
	if (is.character(loc)) {
		if (loc == "right") {
			p <- c(usr[2]+dxy[1], usr[2]+2*dxy[1], usr[3], usr[4])
		} else if (loc == "left") {
			s <- .line.usr(trunc(par("mar")[2]), 2)
			p <- c(s+4*dxy[1], s+5*dxy[1], usr[3], usr[4])
		} else if (loc == "bottom") {
			s <- .line.usr(trunc(par("mar")[1]), 1)
			p <- c(usr[1], usr[2], s+2*dxy[2], s+3*dxy[2])
		} else if (loc == "top") {
			p <- c(usr[1], usr[2], usr[4]+dxy[2], usr[4]+2*dxy[2])
		} else {
			stop(loc)
		}
	}
	x$leg$ext <- p
	x$leg$user <- FALSE
	.get.leg.coords(x)
}

#.plt(r, leg="left", mar=c(3,10,3,3))


.plot.legend <- function(x) {
	if (is.null(x$leg$ext)) {
		x <- .get.leg.extent(x)
	} else {
		x <- .get.leg.coords(x)	
	}
	if (x$leg$type == "continuous") {
		.cont.legend(x)
	}
}
#.plt(r, leg="left", mar=c(3,10,3,3), leg.ext=e)



.plotit <- function(x, leg.ext=NULL, leg.levels=NULL, leg.at=NULL, minmax=NULL, ...) {

	graphics::par(mar=x$mar)	
	
	plot(x$ext[1:2], x$ext[3:4], type = "n", xlab = "", ylab = "", asp=x$asp, ...)
	rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	
	if (x$leg$legend) {	
		if (inherits(leg.ext, "SpatExtent")) {
			leg.ext <- as.data.frame(rbind(as.vector(e)))
		} else if (is.numeric(leg.ext)) {
			leg.ext <- data.frame(rbind(leg.ext))
			if (ncol(leg.ext) != 4) {
				leg.ext <- NULL
			} else {
				colnames(leg.ext) <- c("xmin", "xmax", "ymin", "ymax")
			}
		}
		x$leg$ext <- leg.ext
		if (is.null(leg.levels)) {
			leg.levels <- 5
		} else {
			x$leg$levels <- leg.levels
		}
		x$leg$at <- leg.at
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
	out$leg$shrink <- c(0,0)
	out$leg$main <- NULL
	out$leg$main.cex <-  1	
	
	out
}


.get.plot.data <- function(x, type, cols, maxcell, mar, leg,  ...) {
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
	if (is.na(leg) || isFALSE(leg)) {
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
		x <- .get.plot.data(object, type=type, cols=col, mar=mar, leg=legend, ...)
		x$interpolate <- interpolate

		#y <- .get.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), leg="bottom")
		#y$interpolate <- F

		.plotit(x, ...)
	}
#}

#.plt(r, leg="top", mar=c(2,2,2,2), leg.ext=e, leg.levels=3, leg.at=c(666,999), minmax=c(0,2000))
