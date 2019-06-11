

.getMar <- function(dmar) {
	mar <- graphics::par()$mar
	if (!is.null(.terra_environment$mar)) {
		if (all(mar == (.terra_environment$mar + .terra_environment$dmar))) {
			mar <- .terra_environment$mar
		} 
	}
	.terra_environment$mar <- mar
	.terra_environment$dmar <- dmar
	mar + dmar
}

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


setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function (x, y, cols, maxpixels = 100000, leg.mar, leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1, leg.ext, digits, useRaster = TRUE, zlim, xlab="", ylab="", ...) {

		asp <- 1
		if (couldBeLonLat(x, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} 
		
		if (missing(cols)) cols <- rev(grDevices::terrain.colors(25))
		object <- sampleRegular(x[[y]], maxpixels)
		
		Y <- yFromRow(object, nrow(object):1)
		Z <- t(as.matrix(object, TRUE)[nrow(object):1, , drop = FALSE])
		X <- xFromCol(object, 1:ncol(object))

		if (missing(zlim)) {
			zlim <- range(Z, na.rm=TRUE)
		} else {
			zlim <- sort(zlim)
			Z[Z < zlim[1]] <- zlim[1]
			Z[Z > zlim[2]] <- zlim[2]
		}
#		dv <- dev.list()
#		if (!is.null(dv)) {
#			# excluding pdf, png, ....
#			if (names(dv[length(dv)]) %in% 
#			c("windows", "cairox11", "RStudioGD")){
#				graphics::plot.new()
#			}
#		}

		if (missing(leg.mar)) {
			if (missing(leg.ext)) {
				leg.mar=3
			} else {
				leg.mar=0
			}
		}
		leg.hor <- FALSE
		if (leg.hor) {
			graphics::par(mar=.getMar(c(leg.mar, 0, 0, 0)))
		} else {
			graphics::par(mar=.getMar(c(0, 0, 0, leg.mar)))
		}		
		image(X, Y, Z, col=cols, useRaster=useRaster, asp=asp, xlab=xlab, ylab=ylab, ...)


		if (missing(digits)) {
			dif <- diff(zlim)
			if (dif == 0) {
				digits = 0;
			} else {
				digits <- max(0, -floor(log10(dif/10)))
			}
		}

		u <- unique(na.omit(as.vector(Z)))
		if (length(u) < 10) {
			cols <- .sampleColors(cols, length(u))
		}
		
		usr <- graphics::par()$usr
		dx <- graphics::par()$cxy[1] * graphics::par("cex")	

		if (!missing(leg.ext)) {
			if (inherits(leg.ext, "SpatExtent")) {
				xex <- as.vector(ext(object))
				leg.ext <- .getLegCoords(NULL, xex, leg.shrink, leg.main)
			} else if (is.vector(leg.ext)) {
				leg.ext <- .getLegCoords(NULL, leg.ext, leg.shrink, leg.main)
			} 
			leg.ext.set <- TRUE
		}
		if (missing(leg.ext)) {
			leg.ext.set <- FALSE
			p <- c(usr[2]+dx, usr[2]+2*dx, usr[3], usr[4])
			xex <- as.vector(ext(object))
			leg.ext <- .getLegCoords(p, xex, leg.shrink, leg.main)
		} 


		u <- unique(na.omit(as.vector(Z)))
		if (length(u) < 10) {
			n <- ifelse(leg.ext.set, length(u), 20)
			.fewClassLegend(leg.ext, u, cols, digits, n)
		} else {
			.contLegend(leg.ext, cols, zlim, digits, leg.levels)
		}
		.legMain(leg.main, leg.ext$xmax, leg.ext$ymax, leg.ext$dy, leg.main.cex)
		
		
#		setHook("before.plot.new", 
#			function(...) {
#				m <- graphics::par()$mar
#				m[4] <- mars[4]
#				graphics::par(mar=m)
#				setHook("before.plot.new", NULL, action="replace")
#			}, action="replace")		
	}
)



setMethod("plot", signature(x="SpatRaster", y="missing"), 
	function(x, y, maxpixels=100000, nc, nr, main, maxnl=16, ...)  {
		nl <- min(nlyr(x), maxnl)
		if (nl == 0) {
			stop('SpatRaster has no cell values')
		}
		if (nl==1) {
			plot(x, 1, ...)
			return(invisible(NULL))
		}
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
		
		old.par <- graphics::par(no.readonly = TRUE) 
		on.exit(graphics::par(old.par))
		graphics::par(mfrow=c(nr, nc), mar=c(2, 2, 2, 4))

		maxpixels=maxpixels/(nl/2)
			
		if (missing(main)) {
			main <- names(x)
		} else {
			main <- rep_len(main, nl)
		}		
			
		for (i in 1:nl) {
			image(x[[i]], main=main[i], ...)
		}
	}
)



setMethod("lines", signature(x="SpatRaster"),
function(x, mx=50000, ...) {
	if(prod(dim(x)) > mx) {
		stop('too many lines')
	}
	v <- as.polygons(x)
	lines(v, ...)
}
)
