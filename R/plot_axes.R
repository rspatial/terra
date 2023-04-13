
.plot.axes <- function(x) {

	if (is.null(x$axs$cex.axis)) {
		x$axs$cex.axis = 1
	}
	x$axs$cex.axis <- x$axs$cex.axis * 0.7
	
	if (is.null(x$axs$mgp)) {
		x$axs$mgp = c(2, .25, 0)
	}
	if (is.null(x$axs$tcl)) {
		x$axs$tcl <- -0.25
	}
	if ((!x$clip) & x$draw_grid) {
		x$axs$tck <- 1
		x$axs$mgp <- c(2, .15, 0)
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

	if (x$clip) {
		usr <- x$lim
	} else {
		usr <- graphics::par("usr")
	}
	
	y <- x$axs
	retro <- isTRUE(y$retro) 
	if (retro && (!x$lonlat)) {
		warn("plot", "'retro' labels can only be used with lonlat data") 
		retro <- FALSE
	}
	y$retro <- y$lab <- y$tick <- NULL
	y$line <- NA
	y$outer <- FALSE
	y$line.lab <- NULL
	if (is.null(y$col)) y$col <- "black"
	if (x$clip) {
		lnpts <- crds(as.points(ext(x$lim)))
		lnpts <- rbind(lnpts[4,], lnpts)
	} else {
		lnpts <- crds(as.points(ext(usr)))
		lnpts <- rbind(lnpts[4,], lnpts)
	}
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

			if (x$clip && x$draw_grid && s == 1) {
				clp <- get.clip()
				if (!is.null(clp)) {
					for (i in seq_along(y$at)) {
						lines(rbind(c(y$at[i], clp[3]), c(y$at[i], clp[4])), col=y$col)
					}
				}
			}
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

			if (x$clip && x$draw_grid && s == 2) {
				clp <- get.clip()
				if (!is.null(clp)) {
					for (i in seq_along(y$at)) {
						lines(rbind(c(clp[1], y$at[i]), c(clp[2], y$at[i])), col=y$col)
					}
				}
			}
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
	
	if (x$xlab != "") {
		posx <- usr[1] + diff(usr[1:2])/2
		text(posx, usr[3], x$xlab, pos=1, offset=x$axs$line.lab, cex=x$axs$cex.lab, xpd=TRUE)
	}
	if (x$ylab != "") {
		posy <- usr[3] + diff(usr[3:4])/2
		text(usr[1], posy, x$ylab, pos=2, offset=x$axs$line.lab, srt=90, cex=x$axs$cex.lab, xpd=TRUE)
	}
	
	x
}

