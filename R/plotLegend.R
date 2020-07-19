

#.getMar <- function(dmar) {
#	mar <- graphics::par()$mar
#	if (!is.null(.terra_environment$mar)) {
#		if (all(mar == (.terra_environment$mar + .terra_environment$dmar))) {
#			mar <- .terra_environment$mar
#		} 
#	}
#	.terra_environment$mar <- mar
#	.terra_environment$dmar <- dmar
#	mar + dmar
#}

.legMain <- function(leg.main, xmax, ymax, dy, leg.main.cex) {
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

	if (is.null(p)) {
		xmin <- ext[1]
		xmax <- ext[2]
		ymin <- ext[3]
		ymax <- ext[4]
	} else {
		xmin <- p[1]
		xmax <- p[2]
		ymin <- p[3]
		ymax <- p[4]
		ymin <- max(ymin, ext["ymin"])
		ymax <- min(ymax, ext["ymax"])
	}

	leg.shrink <- rep_len(leg.shrink,2)
	if (!is.null(leg.main)) {
		n <- length(leg.main)		
		leg.shrink[2] <- max(leg.shrink[2], (.05*n)) 
	}

	yd <- ymax - ymin
	ymin <- ymin + yd * leg.shrink[1]
	ymax <- ymax - yd * leg.shrink[2]
    dx <- xmax - xmin
	dy <- ymax - ymin

	data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, dx=dx, dy=dy)
}


.contLegend <- function(e, cols, zlim, digits, leg.levels, ...) {

    nc <- length(cols)
	cols <- rev(cols)
    Y <- seq(e$ymin, e$ymax, length.out=nc+1)
    graphics::rect(e$xmin, Y[-(nc + 1)], e$xmax, Y[-1], col=rev(cols), border=NA, xpd=TRUE)
    graphics::rect(e$xmin, e$ymin, e$xmax, e$ymax, border ="black", xpd=TRUE)
	
	zz <- pretty(zlim, n =(leg.levels+1))	
	zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
    ypos <- e$ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * e$dy
    graphics::segments(e$xmin, ypos, e$xmax+e$dx*0.25, ypos, xpd=TRUE)
    text(e$xmax, ypos, formatC(zz, digits=digits, format = "f"), pos=4, xpd=TRUE, ...)
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



.fewClassLegend <- function(e, u, cols, digits, nspaces, ...) {
	u <- sort(u)
	n <- length(u)
	step <- e$dy / nspaces
    Y <- e$ymax - 0:(n-1) * (step *1.5)
	#to put in the middle
	#mid <- Y[trunc((length(Y)+1)/2)]
	#shift <- max(0, mid - ((e$ymax - e$ymin) / 2))
	#Y <- Y - shift
	
	for (i in 1:n) {
		graphics::rect(e$xmin, Y[i], e$xmax, Y[i]-step, col=cols[i], border="black", xpd=TRUE)
	}
	if (all(u == round(u))) digits = 0
    text(e$xmax, Y-0.5*step, formatC(u, digits=digits, format = "f"), pos=4, xpd=TRUE, ...)
}

.factorLegend <- function(e, u, cols, labs, nspaces, ...) {
	n <- length(u)
	step <- e$dy / nspaces
    Y <- e$ymax - 0:(n-1) * (step *1.5)
	for (i in 1:n) {
		graphics::rect(e$xmin, Y[i], e$xmax, Y[i]-step, col=cols[i], border="black", xpd=TRUE)
	}
    text(e$xmax, Y-0.5*step, labs, pos=4, xpd=TRUE, ...)
}


