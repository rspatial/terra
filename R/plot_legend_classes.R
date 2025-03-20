


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
	bty="n", border="black", seg.len=1, plotlim, yshift=NULL, order=FALSE, sort=FALSE, reverse=FALSE,
	title=NULL, leg_i=1, title.x=NULL, title.y=NULL, title.adj=0.5, title.pos=NULL, 
	title.cex=cex[1], title.col=par("col"), title.font=NULL, ...,
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
	
	if ((!is.null(title.x)) && (!is.null(title.y))) {
		text(x=title.x, y=title.y, labels=title, pos=title.pos, cex=title.cex, xpd=NA, adj=title.adj, font=title.font, col=title.col)
		title <- ""
	}
	
	if (reverse) {
		fill <- rev(fill)
		legend <- rev(legend)
		lty <- rev(lty)
		lwd <- rev(lwd)
		bty <- rev(bty)
		pch <- rev(pch)
		pt.cex <- rev(pt.cex)
		pt.bg <- rev(pt.bg)
		pt.lwd <- rev(pt.lwd)
		cex <- rev(cex)
		density <- rev(density*2)
		angle <- rev(angle)
		border <- rev(border)		
		seg.len <- rev(seg.len)
	}
	
	
	

#points(leg$rect$left+leg$rect$w, leg$rect$top-leg$rect$h, xpd=T)	
	if (grepl("points", geomtype)) {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, plot=FALSE, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font,...)$rect
			xy <- get_legxy(r, plotlim, x, yshift)
			leg <- legend(xy[1], xy[2], legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font,...)
		} else {
			leg <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, pch=pch, pt.cex=pt.cex, pt.bg=pt.bg, pt.lwd=pt.lwd, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)
		}
	} else if (geomtype == "lines") {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, plot=FALSE, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)$rect
			xy <- get_legxy(r, plotlim, x, yshift)
			leg <- legend(xy[1], xy[2], legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)
		} else {
			leg <- legend(x, y, legend, col=fill, xpd=xpd, bty=bty, cex=cex, lty=lty, lwd=lwd, seg.len=seg.len, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)
		}
	} else {
		if (inherits(x, "character")) {
			r <- legend(x, y, legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, plot=FALSE, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)$rect
			xy <- get_legxy(r, plotlim, x, yshift)
			leg <- legend(xy[1], xy[2], legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)
		} else {
			leg <- legend(x, y, legend, fill=fill, xpd=xpd, bty=bty, cex=cex, density=density*2, angle=angle, border=border, title=title, title.adj=title.adj, title.cex=title.cex, title.col=title.col, title.font=title.font, ...)
		}
	}
		
	leg
}

