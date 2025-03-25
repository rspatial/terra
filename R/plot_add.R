
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


add_grid <- function(nx=NULL, ny=nx, col="lightgray", lty="dotted", lwd=1) {

	p <- get.clip()

	## adapted from graphics::grid 
	g.grid.at <- function (side, n, axp, usr2) {
		if (is.null(n)) {
			stopifnot(is.numeric(ax <- axp), length(ax) == 3L)
			graphics::axTicks(side, axp=ax, usr=usr2, log=FALSE)
		}
		else if (!is.na(n) && (n <- as.integer(n)) >= 1L) {
			at <- seq.int(usr2[1L], usr2[2L], length.out = n + 1L)
			at[-c(1L, n + 1L)]
		}
	}

    atx <- if (is.null(nx) || (!is.na(nx) && nx >= 1)) 
        g.grid.at(1L, nx, axp = graphics::par("xaxp"), usr2 = p[1:2])
    aty <- if (is.null(ny) || (!is.na(ny) && ny >= 1)) 
        g.grid.at(2L, ny, axp = graphics::par("yaxp"), usr2 = p[3:4])
    graphics::abline(v = atx, h = aty, col = col, lty = lty, lwd = lwd)
    invisible(list(atx = atx, aty = aty))
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

