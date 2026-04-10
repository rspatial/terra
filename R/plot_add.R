

.txt.loc <- function(x, main=TRUE) {
	if (main) {
		cex <- "cex.main"
		loc <- "loc.main"
	} else {
		cex <- "cex.sub"
		loc <- "loc.sub"
	}
	if (isTRUE(x$clip)) {
		dxy <- graphics::par("cxy") * x[[cex]]
		if (grepl("right", x[[loc]])) {
			px <- x$lim[2]
			pos <- 2
		} else {
			px <- x$lim[1]
			pos <- 4	
		}
		if (grepl("bottom", x[[loc]])) {
			py <- x$lim[3] + dxy[2]/2
		} else {
			py <- x$lim[4] - dxy[2]/2
		}
	} else {
		dxy <- graphics::par("cxy") * x[[cex]]
		usr <- graphics::par("usr")
		if (grepl("right", x[[loc]])) {
			px <- usr[2]
			pos <- 2
		} else {
			px <- usr[1]
			pos <- 4	
		}
		if (grepl("bottom", x[[loc]])) {
			py <- usr[3] + dxy[2]/2
		} else {
			py <- usr[4] - dxy[2]/2
		}
	}
	out <- c(px, py, pos)
	names(out) <- NULL
	out
}


legend_cont <- function(x="right", y=NULL, legend, col,
	size=NULL, title=NULL, at=NULL, digits=NULL, ...) {

	if (inherits(legend, "SpatRaster")) {
		mm <- minmax(legend)
		rng <- c(mm[1,1], mm[2,1])
	} else {
		rng <- range(legend, na.rm=TRUE)
	}
	stopifnot(length(rng) == 2)

	if (missing(col)) col <- .default.pal()
	col <- grDevices::rgb(t(grDevices::col2rgb(col)), maxColorValue=255)

	out <- list()
	out$cols <- col
	out$range <- rng
	out$levels <- 5
	out$fill_range <- TRUE
	out$range_filled <- c(FALSE, FALSE)
	e <- unlist(get.clip())
	usr <- graphics::par("usr")
	tol <- max(abs(usr)) * .Machine$double.eps * 100
	if (!is.null(e) &&
		e[1] >= usr[1] - tol && e[2] <= usr[2] + tol &&
		e[3] >= usr[3] - tol && e[4] <= usr[4] + tol) {
		out$clip <- TRUE
		out$lim <- e[1:4]
	} else {
		out$clip <- FALSE
		out$lim <- usr
	}
	out$axs <- list(lab=NULL)

	out$leg <- list(...)
	horiz <- isTRUE(out$leg$horiz) || (is.character(x) && x %in% c("top", "bottom"))
	if (is.numeric(x)) {
		if (length(x) == 2) {
			out$leg$x <- sort(x)
		} else if (horiz && length(x) == 1) {
			w <- diff(out$lim[1:2]) * 0.6
			out$leg$x <- c(x, x + w)
		} else {
			out$leg$x <- x
		}
	} else {
		out$leg$x <- x
	}
	if (is.numeric(x) && !is.null(y)) {
		if (length(y) == 2) {
			out$leg$y <- sort(y)
		} else if (!horiz && length(y) == 1) {
			h <- diff(out$lim[3:4]) * 0.6
			out$leg$y <- c(y - h, y)
		} else {
			out$leg$y <- y
		}
	} else {
		out$leg$y <- y
	}
	if (!is.null(size)) out$leg$size <- size
	if (!is.null(title)) out$leg$title <- title
	if (!is.null(at)) out$leg$at <- at
	if (is.null(out$leg$reverse)) out$leg$reverse <- FALSE
	if (is.null(out$leg$digits)) {
		if (!is.null(digits)) {
			out$leg$digits <- digits
		} else {
			dif <- diff(rng)
			if ((dif == 0) || (length(dif) == 0)) {
				out$leg$digits <- 0
			} else {
				out$leg$digits <- max(0, -floor(log10(dif / 10)))
			}
		}
	}

	.plot.cont.legend(out)
	invisible(out)
}


add_legend <- function(x, y, xpd=TRUE, ...) {
	if (inherits(x, "character")) {
		e <- unlist(get.clip())
		if (!is.null(e)) {
			rct <- graphics::legend(x=x, y=y, plot=FALSE, xpd=xpd, ...)$rect
			xy <- get_legxy(rct, e[1:4], x, NULL)
			graphics::legend(x=xy[1], y=xy[2], xpd=xpd, ...)
		} else {
			graphics::legend(x=x, y=y, xpd=xpd, ...)
		}
	} else {
		graphics::legend(x=x, y=y, xpd=xpd, ...)
	}
}


add_box <- function(...) {
	e <- unlist(get.clip())
	if (!is.null(e)) {
		bx <- rbind(
			cbind(e[1], e[3:4]),
			cbind(e[2], e[4:3]),
			cbind(e[1], e[3])
		)
		if (is.null(list(...)$xpd)) {
			lines(bx, xpd=TRUE, ...)
		} else {
			lines(bx, ...)		
		}
	}
}

## adapted from graphics::grid 
g.grid.at <- function (side, n, axp, usr2) {
	if (is.null(n)) {
		stopifnot(is.numeric(ax <- axp), length(ax) == 3L)
		graphics::axTicks(side, axp=ax, usr=usr2, log=FALSE)
	} else if (!is.na(n) && (n <- as.integer(n)) >= 1L) {
		at <- seq.int(usr2[1L], usr2[2L], length.out = n + 1L)
		at[-c(1L, n + 1L)]
	}
}

add_grid <- function(nx=NULL, ny=nx, col="lightgray", lty="dotted", lwd=1) {

	p <- unlist(get.clip())
	reset.clip()

    atx <- if (is.null(nx) || (!is.na(nx) && nx >= 1)) 
        g.grid.at(1L, nx, axp=graphics::par("xaxp"), usr2 = p[1:2])
    aty <- if (is.null(ny) || (!is.na(ny) && ny >= 1)) 
        g.grid.at(2L, ny, axp = graphics::par("yaxp"), usr2 = p[3:4])
    graphics::abline(v = atx, h = aty, col = col, lty = lty, lwd = lwd)
    invisible(list(atx = atx, aty = aty))
}

add_abline <- function(h=NULL, v=NULL, ...) {
	p <- unlist(get.clip())
	reset.clip()
	if (!is.null(h)) {
		lines(as.lines(cbind(p[1], p[2], h, h)), ...)
	}
	if (!is.null(v)) {
		lines(as.lines(cbind(v, v, p[3], p[4])), ...)	
	}
}

	
add_mtext <- function(text, side=3, line=0, ...) {

	stopifnot(side %in% 1:4)

	p <- unlist(get.clip())
	h <- graphics::strheight(text, units = "user", ...)	
	
	srt <- 0
	if (side==1) {
		x <- mean(p[1:2])
		y <- p[3] - h - line * h
	} else 	if (side==2) {
		x <- p[1] -1.25 * h - line * h
		y <- mean(p[3:4])
		srt <- 90
	} else 	if (side==3) {
		x <- mean(p[1:2])
		y <- p[4] + h + line * h
	} else {
		x <- p[2] + 1.25 * h + line * h
		y <- mean(p[3:4])
		srt <- 270
	} 
	
	text(x=x, y=y, labels=text, xpd=TRUE, srt=srt, ...)
}


plot_main <- function(x) {
		if (isTRUE(x$main != "")) {
			pos <- 3
			if (is.null(x$loc.main)) {
				if (isTRUE(x$clip)) {
					x$loc.main <- c(x$lim[1] + diff(x$lim[1:2]) / 2, x$lim[4])
				} else {
					usr <- graphics::par("usr")			
					x$loc.main <- c(usr[1] + diff(usr[1:2]) / 2, usr[4])			
				}
			} else if (inherits(x$loc.main, "character")) {
				xyp <- .txt.loc(x)
				x$loc.main <- xyp[1:2]
				pos <- xyp[3]
			}	
			if (isTRUE(x$halo.main)) {
				.halo(x$loc.main[1], x$loc.main[2], x$main, pos=pos, offset=x$line.main, cex=x$cex.main, 
					font=x$font.main, col=x$col.main, xpd=TRUE, hc=x$halo.main.hc, hw=x$halo.main.hw)
			} else {
				text(x$loc.main[1], x$loc.main[2], x$main, pos=pos, offset=x$line.main, cex=x$cex.main, 
					font=x$font.main, col=x$col.main, xpd=TRUE)
			}
		}
		if (isTRUE(x$sub != "")) { 
			pos <- 1
			if (is.null(x$loc.sub)) {
				if (isTRUE(x$clip)) {
					x$loc.sub <- c(x$lim[1] + diff(x$lim[1:2]) / 2, x$lim[3])
				} else {
					usr <- graphics::par("usr")			
					x$loc.sub <- c(usr[1] + diff(usr[1:2]) / 2, usr[3])			
				}
			} else if (inherits(x$loc.sub, "character")) {
				xyp <- .txt.loc(x, FALSE)
				x$loc.sub <- xyp[1:2]
				pos <- xyp[3]
			}
			if (isTRUE(x$halo.main)) {
				.halo(x$loc.sub[1], x$loc.sub[2], x$sub, pos=pos, offset=x$line.sub, cex=x$cex.sub, 
					font=x$font.sub, col=x$col.sub, xpd=TRUE, hc=x$halo.sub.hc, hw=x$halo.sub.hw)
			} else {
				text(x$loc.sub[1], x$loc.sub[2], x$sub, pos=pos, offset=x$line.sub, cex=x$cex.sub, 
					font=x$font.sub, col=x$col.sub, xpd=TRUE)
			}
		}
}

