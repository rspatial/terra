

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

	if (x$leg$horizontal) {
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
			p <- c(xmin, xmax, s+1.75*dxy[2], s+2.25*dxy[2])
		} else if (loc == "top") {
			p <- c(xmin, xmax, ymax+.6*dxy[2], ymax+1*dxy[2])
		} else if (loc == "topright") {
			p <- c(xmax+dxy[1], xmax+2*dxy[1], ymin + (ymax-ymin) / 2, ymax)
			x$leg$x <- "right"
		} else if (loc == "bottomright") {
			p <- c(xmax+dxy[1], xmax+2*dxy[1], ymin, ymax - (ymax-ymin) / 2)
			x$leg$x <- "right"
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




.plot.cont.legend <- function(x, ...) {

	if (!is.null(x$leg[["tic"]])) {
		accepted <- c("in", "out", "none", "middle", "through", "throughout")
		tics <- accepted[pmatch(x$leg$tic[1], accepted[-6], 6)]
	} else {
		tics <- "throughout"
	}

	if (!is.null(x$leg[["tic.box.col"]])) {
		ticboxcol <- x$leg$tic.box.col[1]
	} else {
		ticboxcol <- "black"
	}
	if (!is.null(x$leg[["tic.col"]])) {
		ticcol <- x$leg$tic.col[1]
	} else {
		ticcol <- "black"
	}
	if (!is.null(x$leg[["tic.lwd"]])) {
		ticlwd <- x$leg$tic.lwd[1]
	} else {
		ticlwd <- 1
	}
	boxlwd <- 1 # lwdd?

	
	x$leg$horizontal <- FALSE
	if (is.null(x$leg[["x"]])) {
		x$leg$x <- "right"
	} else if (!(x$leg$x %in% c("left", "right", "top", "bottom", "topright", "bottomright"))) {
		x$leg$x <- "right"	
	} else if (x$leg$x %in% c("top", "bottom")) {
		x$leg$horizontal <- TRUE
	}

	if (is.null(x$leg[["ext"]])) {
		x <- .get.leg.extent(x)
	} else {
		x <- .get.leg.coords(x)
	}

	if (is.null(x$leg[["cex"]])) {
		cex <- 1
	} else {
		cex <- x$leg$cex
	}
	cex <- cex * 0.8
	
	if (!is.null(x$leg[["rotate"]])) {
		srt <- ifelse(isTRUE(x$leg$rotate), 90, 0)
	} else if (!is.null(x$leg$srt)) {
		srt <- x$leg$srt
	} else {
		srt <- 0	
	}
	
	cols <- rev(x[["cols"]])
	nc <- length(cols)

	zlim <- x[["range"]]
	zz <- x$leg[["at"]]
	if (is.null(zz)) {
		if (is.null(x[["levels"]])){
			x$levels <- 5
		}
		zz <- pretty(zlim, n =(x$levels+1))
		zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
	}
	zztxt <- x$leg[["labels"]]
	if (is.null(zztxt)) {
		zztxt <- formatC(zz, digits=x$leg[["digits"]], format = "f")
		if (x$fill_range) {
			if (isTRUE(x$range_filled[1])) zztxt[1] <- paste0("< ", zztxt[1])		
			if (isTRUE(x$range_filled[2])) zztxt[length(zztxt)] <- paste0("> ", zztxt[length(zztxt)])		
		}
	}
	e <- x$leg[["ext"]]
	if (x$leg$x %in% c("left", "right")) {
		Y <- seq(e$ymin, e$ymax, length.out=nc+1)
		graphics::rect(e$xmin, Y[-(nc + 1)], e$xmax, Y[-1], col=rev(cols), border=NA, xpd=NA, lwd=boxlwd)
		ypos <- e$ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dy
		if (x$leg$x == "right") {
			if (tics == "throughout") {
				graphics::segments(e$xmin, ypos, e$xmax+e$dx*0.25, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "through") {
				graphics::segments(e$xmin, ypos, e$xmax, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "middle") {
				graphics::segments(e$xmin+e$dx*0.25, ypos, e$xmax-e$dx*0.25, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "in") {
				graphics::segments(e$xmax-e$dx*0.25, ypos, e$xmax, ypos, xpd=NA, col=ticcol, lwd=ticlwd)			
			} else if (tics == "out") {
				graphics::segments(e$xmax, ypos, e$xmax+e$dx*0.25, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			}
			text(e$xmax, ypos, zztxt, pos=4, xpd=NA, cex=cex, srt=srt, ...)
		} else {
			if (tics == "throughout") {
				graphics::segments(e$xmin-e$dx*0.25, ypos, e$xmax, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "through") {
				graphics::segments(e$xmin, ypos, e$xmax, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "middle") {
				graphics::segments(e$xmin+e$dx*0.25, ypos, e$xmax-e$dx*0.25, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "in") {
				graphics::segments(e$xmin, ypos, e$xmin+e$dx*0.25, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "out") {
				graphics::segments(e$xmin-e$dx*0.25, ypos, e$xmin, ypos, xpd=NA, col=ticcol, lwd=ticlwd)
			}
			text(e$xmin, ypos, zztxt, pos=2, xpd=NA, cex=cex, srt=srt, ...)
		}
	} else { # top, bottom
		X <- seq(e$xmin, e$xmax, length.out=nc+1)
		graphics::rect(X[-(nc + 1)], e$ymin, X[-1], e$ymax, col=rev(cols), border=NA, xpd=NA, lwd=boxlwd)
		xpos <- e$xmin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dx
		
		if (x$leg$x == "bottom") {
			if (tics == "throughout") {
				graphics::segments(xpos, e$ymin-e$dy*0.25, xpos, e$ymax, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "through") {
				graphics::segments(xpos, e$ymin, xpos, e$ymax, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "middle") {
				graphics::segments(xpos, e$ymax-e$dy*0.25, xpos, e$ymin+e$dy*0.25, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "in") {
				graphics::segments(xpos, e$ymin+e$dy*0.25, xpos, e$ymin, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "out") {
				graphics::segments(xpos, e$ymin-e$dy*0.25, xpos, e$ymin, xpd=NA, col=ticcol, lwd=ticlwd)
			}
			text(xpos, e$ymin-e$dy, zztxt, pos=NULL, xpd=NA, cex=cex, srt=srt, ...)
		} else {
			if (tics == "throughout") {
				graphics::segments(xpos, e$ymin, xpos, e$ymax+e$dy*0.25, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "through") {
				graphics::segments(xpos, e$ymin, xpos, e$ymax, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "middle") {
				graphics::segments(xpos, e$ymax-e$dy*0.25, xpos, e$ymin+e$dy*0.25, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "in") {
				graphics::segments(xpos, e$ymax, xpos, e$ymax-e$dy*0.25, xpd=NA, col=ticcol, lwd=ticlwd)
			} else if (tics == "out") {
				graphics::segments(xpos, e$ymax, xpos, e$ymax+e$dy*0.25, xpd=NA, col=ticcol, lwd=ticlwd)
			}
			text(xpos, e$ymax+1.5*e$dy, zztxt, pos=NULL, xpd=NA, cex=cex, srt=srt, ...)
		}
	}
	graphics::rect(e$xmin, e$ymin, e$xmax, e$ymax, border=ticboxcol, xpd=NA)


    if (isTRUE("title" %in% names(x$leg))) {
		leg_i <- x$leg$leg_i
		if (is.null(leg_i)) leg_i = 1
	    if (leg_i <= length(x$leg$title)) {
			legtitle <- x$leg$title[leg_i]
		} else {
			legtitle <- x$leg$title[1]		
		}
		if ((is.null(x$leg[["title.x"]])) || (is.null(x$leg[["title.y"]]))) {
			e <- x$leg$ext
			x$leg$title.y <- e$ymax
			if (x$leg$horizontal) {
				x$leg$title.x <- e$xmin + (e$xmax - e$xmin) / 2
				if (x$leg$x	== "top") {
					x$leg$title.y <- e$ymax + 2 * (e$ymax - e$ymin)
				}
			} else {
				x$leg$title.x <- e$xmax	
			}
			if (length(legtitle) > 1) { # or perhaps !inherits(legtitle, "expression")
				if (x$leg$horizontal) {
					legtitle <- paste(legtitle, collapse=" ")
				} else {
					legtitle <- paste(legtitle, collapse="\n")		
				}
			} 
		}
		pos <- 3
		if (!is.null(x$leg[["title.pos"]])) pos <- x$leg$title.pos

		# offset=.5*graphics::strheight("a",cex=x$leg$title.cex)
		text(x=x$leg$title.x, y=x$leg$title.y, labels=legtitle, pos=pos, cex=x$leg$title.cex, xpd=NA, adj=x$leg$title.adj, font=x$leg$title.font, col=x$leg$title.col, srt=x$leg$title.srt)
	}
	x
}


