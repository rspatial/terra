
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



.get.leg.coords <- function(x) {

	if (is.null(x$leg$ext)) {
		if (x$clip) {
			p <- x$leg$ext <- x$lim
		} else {
			p <- x$leg$ext <- graphics::par("usr")		
		}
	} else {
		p <- as.vector(x$leg$ext)
	}
	xmin <- p[1]
	xmax <- p[2]
	ymin <- p[3]
	ymax <- p[4]
	flip <- FALSE
	
	if (!is.null(x$leg$shrink)) {
		s <- x$leg$shrink
		if ((s[1] <= 1) & (s[1] >= 0.5)) {
			s[1] <- 2*(s[1] - 0.5)
		} else if (s[1] < 0.5) {
			s[1] <- (2*(0.5 - s[1]))
			flip <- TRUE
		}
		x$leg$size <- s
	} 
	
	if (is.null(x$leg$size)) {
		x$leg$size <- c(1,1)
	} else if (length(x$leg$size) == 1) {
		x$leg$size <- c(x$leg$size, 1)
	}
	if (x$leg$size[1] < 0) flip <- TRUE
	x$leg$size <- abs(x$leg$size)

	if (!is.null(x$leg$main)) {
		n <- length(x$leg$main)
		x$leg$size[1] <- min(x$leg$size[1], (1 - .05*n))
	}

	horiz <- isTRUE(x$leg$x %in% c("top", "bottom"))
	if (horiz) {
#		xd <- (xmax - xmin) * x$leg$size[2]
#		xmin <- xmin + xd 
#		xmax <- xmax - xd

		rhalf <- (xmax - xmin) / 2
		xmid <- xmin + rhalf
		xd <- rhalf * x$leg$size[1]
		xmin <- xmid - xd 
		xmax <- xmid + xd
		
#		yd <- (ymax - ymin) * x$leg$size[1]/1.5
#		ymin <- ymin + yd
#		ymax <- ymax - yd

		yd <- ymax - ymin
		if (x$leg$x == "top") {
			ymax <- ymin + yd * x$leg$size[2] 		
		} else {
			ymin <- ymax - yd * x$leg$size[2] 
		}
		if (flip) {
			tmp <- xmin
			xmin <- xmax
			xmax <- tmp
		}
	} else {

		rhalf <- (ymax - ymin) / 2
		ymid <- ymin + rhalf
		yd <- rhalf * x$leg$size[1]
		ymin <- ymid - yd 
		ymax <- ymid + yd

		xd <- xmax - xmin
		#xmin <- xmin + xd * x$leg$size[2]/5
		#xmax <- xmax - xd * x$leg$size[2]/5
		xmax <- xmin + xd * x$leg$size[2] 
		
		if (flip) {
			tmp <- ymin
			ymin <- ymax
			ymax <- tmp
		}
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
	if (x$clip) {
		usr <- x$lim
	} else {
		usr <- graphics::par("usr")
	}
	xmin <- usr[1]
	xmax <- usr[2]
	ymin <- usr[3]
	ymax <- usr[4]
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
			p <- c(xmin, xmax, s+1.75*dxy[2], s+2.5*dxy[2])
		} else if (loc == "top") {
			p <- c(xmin, xmax, ymax+dxy[2], ymax+1.75*dxy[2])
		} else { #if (loc == "right" or "default" 
			p <- c(xmax+dxy[1], xmax+2*dxy[1], ymin, ymax)
			if (isTRUE(x$leg$yshift)) {
				hy <- (ymax - ymin) / 2
				p[3:4] <- p[3:4] - hy
			}
		}
	}
	x$leg$ext <- p
	x$leg$user <- FALSE
	.get.leg.coords(x)
}

.txt.loc <- function(x) {
	if (isTRUE(x$clip)) {
		dxy <- graphics::par("cxy") * x$cex.main
		if (grepl("right", x$loc.main)) {
			px <- x$lim[2]
			pos <- 2
		} else {
			px <- x$lim[1]
			pos <- 4	
		}
		if (grepl("bottom", x$loc.main)) {
			py <- x$lim[3] + dxy[2]/2
		} else {
			py <- x$lim[4] - dxy[2]/2
		}
	} else {
		dxy <- graphics::par("cxy") * x$cex.main
		usr <- graphics::par("usr")
		if (grepl("right", x$loc.main)) {
			px <- usr[2]
			pos <- 2
		} else {
			px <- usr[1]
			pos <- 4	
		}
		if (grepl("bottom", x$loc.main)) {
			py <- usr[3] + dxy[2]/2
		} else {
			py <- usr[4] - dxy[2]/2
		}
	}
	out <- c(px, py, pos)
	names(out) <- NULL
	out
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
	if (is.null(cex)) cex <- 1
	cex <- cex * 0.8
	
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


    if (!is.null(x$leg$title)) {
		leg_i <- x$leg$leg_i
		if (is.null(leg_i)) leg_i = 1
	    if (leg_i <= length(x$leg$title)) {
			legtitle <- x$leg$title[leg_i]
		} else {
			legtitle <- x$leg$title[1]		
		}
		e <- x$leg$ext
		if (x$leg$x %in% c("top", "bottom")) {
			txt <- paste(legtitle, collapse=" ")
		} else {
			txt <- paste(legtitle, collapse="\n")		
		}
		# offset=.5*graphics::strheight("a",cex=x$leg$title.cex)
		text(x=e$xmax, y=e$ymax, labels=txt, pos=3, cex=x$leg$title.cex, xpd=NA)
	}
	x
}

get_legxy <- function(r, e, pos, yshift) {
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
	
	if (!is.null(yshift)) {
		hy <- (e[4] - e[3]) / 2
		xy[2] <- xy[2] - hy
	}
	xy
}


.plot.class.legend <- function(x, y, legend, fill, xpd=NA, cex=1, geomtype="",
	lty=1, lwd=1, pch=1, angle=45, density=NULL, pt.cex = 1, pt.bg="black", pt.lwd=1, 
	bty="n", border="black", seg.len=1, plotlim, yshift=NULL, title=NULL, leg_i=1, ...,
# catch and kill
	merge, trace, size) {

	cex <- cex * 0.8
	if (x %in% c("top", "default")) {
		#usr <- graphics::par("usr")
		x <- plotlim[2]
		y <- plotlim[4]
	}
	
	if (is.null(leg_i)) leg_i = 1
    if (leg_i <= length(title)) {
		title <- title[leg_i]
	} else {
		title <- title[1]		
	}
#points(leg$rect$left+leg$rect$w, leg$rect$top-leg$rect$h, xpd=T)	
	if (grepl("points", geomtype)) {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, plot=FALSE, title=title, ...)$rect
			xy <- get_legxy(r, plotlim, x, yshift)
			leg <- legend(xy[1], xy[2], legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, title=title, ...)
		} else {
			leg <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, title=title,...)
		}
	} else if (geomtype == "lines") {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, plot=FALSE, title=title,, ...)$rect
			xy <- get_legxy(r, plotlim, x, yshift)
			leg <- legend(xy[1], xy[2], legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, title=title,, ...)
		} else {
			leg <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, title=title, ...)
		}
	} else {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, plot=FALSE, title=title, ...)$rect
			xy <- get_legxy(r, plotlim, x, yshift)
			leg <- legend(xy[1], xy[2], legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, title=title, ...)
		} else {
			leg <- legend(x, y, legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, title=title, ...)
		}
	}
	leg
}


add_legend <- function(x, y, ...) {
	if (inherits(x, "character")) {
		e <- unlist(get.clip())
		if (!is.null(e)) {
			rct <- graphics::legend(x=x, y=y, plot=FALSE, ...)$rect
			xy <- get_legxy(rct, e[1:4], x, NULL)
			graphics::legend(x=xy[1], y=xy[2], ...)
		} else {
			graphics::legend(x=x, y=y, ...)
		}
	} else {
		graphics::legend(x=x, y=y, ...)
	}
}


