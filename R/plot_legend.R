
.get_breaks <- function(x, n, method, r=NULL) {
	#x <- x[!is.na(x)]
	if (is.function(method)) {
		if (!is.null(r)) {
			x[(x<r[1]) | (x>r[2])] <- NA
		}
		breaks <- method(x)
	} else if (method[1]=="cases") {
		if (!is.null(r)) {
			x[(x<r[1]) | (x>r[2])] <- NA
		}
		n <- n+1
		i <- seq(0, 1, length.out=n)
		breaks <- quantile(x, i, na.rm=TRUE)
		breaks <- unique(breaks)
		if ((breaks[1] %% 1) != 0) {
			breaks[1] <- breaks[1] - 0.000001
		}
		if ((breaks[n] %% 1) != 0) {
			breaks[n] <- breaks[n] + 0.000001
		}
	} else { # if (method=="eqint") {
		if (is.null(r)) {
			r <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
		}
		small <- 1e-16
		if ((r[1] %% 1) != 0) { r[1] <- r[1] - small }
		if ((r[2] %% 1) != 0) { r[2] <- r[2] + small }
		breaks <- seq(r[1] , r[2], length.out=n+1)
	}
	breaks
}

.get_nrnc <- function(nr, nc, nl) {
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
	c(nr, nc)
}


retro_labels <- function(x, lat=TRUE) {
	if ((is.null(x)) || (!is.numeric(x))) {
		return(x)
	}
	if ((length(x) > 1) && (min(diff(x)) <= 1/120)) {
		d <- floor(x)
		m <- floor(60*(x - d))
		s <- round(3600*(x - d - m/60))
	} else {
		d <- floor(x)
		m <- round(60*(x - d))
		s <- 0
	}

	if (lat) {
		h <- c("S", "", "N")[sign(d)+2]
	} else {
		h <- c("W", "", "E")[sign(d)+2]
	}
	d <- abs(d)
	i <- (s == 0) & (m == 0)
	j <- (s == 0) & (m != 0)

	m <- formatC(m, width=2, flag="0")
	s <- formatC(s, width=2, flag="0")
	r <- paste0(d, "\u00B0" , m, "'", s, '"', h)
	r[i] <- paste0(d[i], "\u00B0" , h[i])
	r[j] <- paste0(d[j], "\u00B0" , m[j], "'", h[j])	
	r
}



.plot.axes <- function(x) {

	if (is.null(x$axs$cex.axis)) {
		x$axs$cex.axis = 0.7
	}
	if (is.null(x$axs$mgp)) {
		x$axs$mgp = c(2, .25, 0)
	}
	if (is.null(x$axs$tcl)) {
		x$axs$tcl <- -0.25
	}
	if (x$draw_grid) {
		x$axs$tck <- 1
		x$axs$mgp = c(2, .15, 0)
	}

	xlab <- ylab <- NULL
	if (!is.null(x$axs$labels)) {
		xlab <- ylab <- x$axs$labels
	}
	if (!is.null(x$axs$xlabs)) {
		xlab <- x$axs$xlabs
		x$axs$xlabs <- NULL
	}
	if (!is.null(x$axs$ylabs)) {
		ylab <- x$axs$ylabs
		x$axs$ylabs <- NULL
	}

	xat <- yat <- NULL
	if (!is.null(x$axs$at)) {
		xat <- yat <- x$axs$at
	}
	if (!is.null(x$axs$xat)) {
		xat <- x$axs$xat
		x$axs$xat <- NULL
	}
	if (!is.null(x$axs$yat)) {
		yat <- x$axs$yat
		x$axs$yat <- NULL
	}

	sides <- unique(x$axs$side)
	if (!is.null(sides)) sides <- round(sides)
	sides <- sides[sides > 0 & sides < 5]
	if (is.null(sides)) {
		x$axs$side <- sides <- 1:2
	}

	ticks <- x$axs$tick 
	if (is.null(ticks)) {
		x$axs$tick <- ticks <- sides
	}
	labs <- x$axs$lab
	if (is.null(labs)) {
		x$axs$lab <- labs <- sides
	} 

#	usr <- graphics::par("usr")
	usr <- x$lim
	y <- x$axs
	retro <- isTRUE(y$retro)
	y$retro <- y$lab <- y$tick <- NULL
	y$line <- NA
	y$outer <- FALSE
	if (is.null(y$col)) y$col <- "black"
	lnpts <- crds(as.points(ext(x$lim)))
	lnpts <- rbind(lnpts[4,], lnpts)
	
	for (s in 1:4) {
		y$side <- s
		y$labels <- NULL
		if (s %in% c(1,3)) {
			ur <- usr[2] - usr[1]
			edg <- c(usr[1]-10*ur, usr[2]+10*ur)
			if (is.null(xat)) {
				axt <- graphics::axTicks(s)
				y$at <- axt[axt >= usr[1] & axt <= usr[2]]
			} else {
				y$at <- xat
			}
			if (is.null(xlab)) {
				y$labels <- if (retro) retro_labels(y$at, lat=FALSE) else y$at
			} else {
				y$labels <- xlab
			}
			y$pos <- ifelse(s==1, usr[3], usr[4])

		} else {
			ur <- usr[4] - usr[3]
			edg <- c(usr[3]-10*ur, usr[4]+10*ur)
			if (is.null(yat)) {
				axt <- graphics::axTicks(s)
				y$at <- axt[axt >= usr[3] & axt <= usr[4]]
			} else {
				y$at <- yat
			}
			if (is.null(ylab)) {
				y$labels <- if (retro) retro_labels(y$at, lat=TRUE) else y$at
			} else {
				y$labels <- ylab
			}
			y$pos <- ifelse(s==2, usr[1], usr[2])
		}
		z <- y
		z$lwd <- 0

		if (s %in% labs) {			
			z$lwd.ticks <- 0
			do.call(graphics::axis, z)
		}
		z$labels <- FALSE
		if (s %in% ticks) {
			z$lwd <- 0
			z$lwd.ticks <- y$lwd.ticks
			if (is.null(z$lwd.ticks)) z$lwd.ticks <- 1
			do.call(graphics::axis, z)
		}
		if (s %in% sides) {
			lin <- lnpts[s:(s+1), ]
			if (is.null(y$lty)) {
				lty <- 1
			} else {
				lty <- y$lty
			}
			lines(lin, y$lwd, lty=lty, col=y$col)
			#d <- diff(edg) * 10
			#z$at <- edg + c(-d, d)
			#z$lwd.ticks <- 0
			#z$lwd <- y$lwd
			#do.call(graphics::axis, z)
		} 
	}
	x
}


.get.leg.coords <- function(x) {

	if (is.null(x$leg$ext)) {
		p <- x$leg$ext <- x$lim
	} else {
		p <- as.vector(x$leg$ext)
	}
	xmin <- p[1]
	xmax <- p[2]
	ymin <- p[3]
	ymax <- p[4]

	if (is.null(x$leg$shrink)) {
		leg.shrink <- c(0,0)
	} else {
		leg.shrink <- rep_len(x$leg$shrink,2)
	}
	if (!is.null(x$leg$main)) {
		n <- length(x$leg$main)
		leg.shrink[2] <- max(x$leg$shrink[2], (.05*n))
	}

	if (isTRUE(x$leg$x=="bottom")) {
		xd <- xmax - xmin
		xmin <- xmin + xd * leg.shrink[1]
		xmax <- xmax - xd * leg.shrink[2]
		yd <- ymax - ymin
		ymin <- ymin + yd * leg.shrink[1]/1.5
		ymax <- ymax - yd * leg.shrink[2]/1.5
	} else {
		yd <- ymax - ymin
		ymin <- ymin + yd * leg.shrink[1]
		ymax <- ymax - yd * leg.shrink[1]
		xd <- xmax - xmin
		xmin <- xmin + xd * leg.shrink[2]/5
		xmax <- xmax - xd * leg.shrink[2]/5
    }
	dx <- xmax - xmin
	dy <- ymax - ymin

	x$leg$ext <- data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, dx=dx, dy=dy)
	x
}


.line.usr <- function(line, side) {
## https://stackoverflow.com/questions/30765866/get-margin-line-locations-in-log-space/30835971#30835971

	lh <- graphics::par("cin")[2] * graphics::par("cex") * graphics::par("lheight")
	x_off <- diff(graphics::grconvertX(c(0, lh), "inches", "npc"))
	y_off <- diff(graphics::grconvertY(c(0, lh), "inches", "npc"))
	if (side == 1) {
		graphics::grconvertY(-line * y_off, "npc", "user")
	} else if (side ==2) {
		graphics::grconvertX(-line * x_off, "npc", "user")
	} else if (side ==3) {
		graphics::grconvertY(1 + line * y_off, "npc", "user")
	} else {
		graphics::grconvertX(1 + line * x_off, "npc", "user")
	}
}

.get.leg.extent <- function(x) {
	#usr <- graphics::par("usr")
	dxy <- graphics::par("cxy") * graphics::par("cex")
	loc <- x$leg$x
	xmin <- x$ext[1]
	xmax <- x$ext[2]
	ymin <- x$ext[3]
	ymax <- x$ext[4]
	p <- NULL
	if (is.character(loc)) {
		if (loc == "left") {
			#s <- .line.usr(trunc(graphics::par("mar")[2]), 2)
			#p <- c(s+4*dxy[1], s+5*dxy[1], ymin, ymax)	
			if (any(2 %in% x$axs$lab)) {
				p <- c(xmin-4*dxy[1], xmin-3*dxy[1], ymin, ymax)			
			} else {
				p <- c(xmin-2*dxy[1], xmin-dxy[1], ymin, ymax)
			}
		} else if (loc == "bottom") {
			s <- .line.usr(trunc(graphics::par("mar")[1]), 1)
			p <- c(xmin, xmax, s+2*dxy[2], s+3*dxy[2])
		} else if (loc == "top") {
			p <- c(xmin, xmax, ymax+dxy[2], ymax+2*dxy[2])
		} else { #if (loc == "right" or "default" 
			p <- c(xmax+dxy[1], xmax+2*dxy[1], ymin, ymax)
		}
	}
	x$leg$ext <- p
	x$leg$user <- FALSE
	.get.leg.coords(x)
}





.leg.main <- function(x) {
	leg <- x$leg
    if (!is.null(leg$title)) {
		e <- leg$ext
		n <- length(leg$title)
		ymax <- e$ymax + 0.05 * e$dy

		for (i in 1:n) {
			if (x$leg$x == "right") {
				text(x=e$xmax, y=ymax+(n-i)*0.05* e$dy,
					labels = leg$title[i], cex = leg$title.cex, xpd=TRUE)
			} else if (x$leg$x == "left") {
				text(x=e$xmin, y=ymax+(n-i)*0.05* e$dy,
					labels = leg$title[i], cex = leg$title.cex, xpd=TRUE)
			} else { # default
				ymax <- e$ymax + e$dy
				text(x=(e$xmin+e$xmax)/2, y=ymax+(n-i)*0.05* e$dy,
					labels = leg$title[i], cex = leg$title.cex, xpd=TRUE)
			}
		}
	}
	x
}


.plot.cont.legend <- function(x, ...) {

	if (is.null(x$leg$x)) {
		x$leg$x <- "right"
	} else if (!(x$leg$x %in% c("left", "right", "top", "bottom"))) {
		x$leg$x <- "right"	
	}

	if (is.null(x$leg$ext)) {
		x <- .get.leg.extent(x)
	} else {
		x <- .get.leg.coords(x)
	}

	cex <- x$leg$cex
	if (is.null(cex)) cex <- 0.8
	rotate <- isTRUE(x$leg$rotate)
	srt <- ifelse(rotate, 90, 0)

	cols <- rev(x$cols)
	nc <- length(cols)

	zlim <- x$range
	zz <- x$leg$at
	if (is.null(zz)) {
		if (is.null(x$levels)){
			x$levels <- 5
		}
		zz <- pretty(zlim, n =(x$levels+1))
		zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
	}
	zztxt <- x$leg$labels
	if (is.null(zztxt)) {
		zztxt <- formatC(zz, digits=x$leg$digits, format = "f")
	}
	e <- x$leg$ext
	if (x$leg$x %in% c("left", "right")) {
		Y <- seq(e$ymin, e$ymax, length.out=nc+1)
		graphics::rect(e$xmin, Y[-(nc + 1)], e$xmax, Y[-1], col=rev(cols), border=NA, xpd=NA)
		ypos <- e$ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dy
		if (x$leg$x == "right") {
			graphics::segments(e$xmin, ypos, e$xmax+e$dx*0.25, ypos, xpd=NA)
			text(e$xmax, ypos, zztxt, pos=4, xpd=NA, cex=cex, ...)
		} else {
			graphics::segments(e$xmin-e$dx*0.25, ypos, e$xmax, ypos, xpd=NA)
			text(e$xmin, ypos, zztxt, pos=2, xpd=NA, cex=cex, ...)
		}
	} else {
		X <- seq(e$xmin, e$xmax, length.out=nc+1)
		graphics::rect(X[-(nc + 1)], e$ymin, X[-1], e$ymax, col=rev(cols), border=NA, xpd=NA)
		xpos <- e$xmin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dx
		if (x$leg$x == "bottom") {
			graphics::segments(xpos, e$ymin-e$dy*0.25, xpos, e$ymax, xpd=NA)
			text(xpos, e$ymin, zztxt, pos=1, xpd=NA, cex=cex)
		} else {
			graphics::segments(xpos, e$ymin, xpos, e$ymax+e$dy*0.25, xpd=NA)
			text(xpos, e$ymax+e$dy*0.25, zztxt, pos=3, xpd=NA, cex=cex)
		}
	}
	graphics::rect(e$xmin, e$ymin, e$xmax, e$ymax, border ="black", xpd=NA)

	x$leg.main <- .leg.main(x)
	x
}

get_legxy <- function(r, e, pos) {
	xy <- c(r$left, r$top)
	if (grepl("top", pos)) {
		xy[2] <- e[4]
	} else if (grepl("bottom", pos)) {
		xy[2] <- e[3] + r$h
	}

	if (grepl("left", pos)) {
		xy[1] <- e[1]
	} else if (grepl("right", pos)) {
		xy[1] <- e[2] - r$w
	}
	xy
}


.plot.class.legend <- function(x, y, legend, fill, xpd=TRUE, cex=0.8, geomtype="",
	lty=1, lwd=1, pch=1, angle=45, density=NULL, pt.cex = 1, pt.bg="black", pt.lwd=1, 
	bty="n", border="black", seg.len=1, plotlim,
# catching
	merge, trace,...) {

	if (x %in% c("top", "default")) {
		#usr <- graphics::par("usr")
		x <- plotlim[2]
		y <- plotlim[4]
	}
	
#points(leg$rect$left+leg$rect$w, leg$rect$top-leg$rect$h, xpd=T)	
	if (grepl("points", geomtype)) {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, plot=FALSE, ...)$rect
			xy <- get_legxy(r, plotlim, x)
			leg <- legend(xy[1], xy[2], legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, ...)
		} else {
			leg <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, ...)
		}
	} else if (geomtype == "lines") {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, plot=FALSE, ...)$rect
			xy <- get_legxy(r, plotlim, x)
			leg <- legend(xy[1], xy[2], legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, ...)
		} else {
			leg <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, ...)
		}
	} else {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, plot=FALSE, ...)$rect
			xy <- get_legxy(r, plotlim, x)
			leg <- legend(xy[1], xy[2], legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, ...)
		} else {
			leg <- legend(x, y, legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, ...)
		}
	}
	leg
}

