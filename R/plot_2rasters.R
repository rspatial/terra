# Author: Robert J. Hijmans
# Date :  June 2019
# Version 1.0
# License GPL v3


.scatterPlotRaster <- function(x, y, maxcell=100000, warn=TRUE, cex, xlab, ylab, nc, nr, 
	maxnl=16, main, add=FALSE, smooth=FALSE, gridded=FALSE, ncol=25, nrow=25, asp=NA,
	colramp=grDevices::colorRampPalette(c("white", grDevices::blues9)), ...) {

	compareGeom(x, y, lyrs=FALSE, crs=FALSE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE)
	nlx <- nlyr(x)
	nly <- nlyr(y)

	maxnl <- max(1, round(maxnl))
	nl <- min(max(nlx, nly), maxnl)
	if (nl > maxnl) {
		nl <- maxnl
		if (nlx > maxnl) {
			x <- x[[1:maxnl]]
			nlx <- maxnl
		}
		if (nly > maxnl) {
			y <- y[[1:maxnl]]
			nly <- maxnl
		}
	}

	if (nlx < nly) {
		x <- x[[rep_len(1:nlx, nly)]]
		nlx <- nly
	} else if (nly < nlx) {
		y <- y[[rep_len(1:nly, nlx)]]	
		nly <- nlx
	}


	if (missing(main)) {
		main <- ""
	}

	if (missing(xlab)) {
		ln1 <- names(x)
	} else {
		ln1 <- xlab
		if (length(ln1) == 1) {
			ln1 <- rep(ln1, nlx)
		}
	}
	if (missing(ylab)) {
		ln2 <- names(y)
	} else {
		ln2 <- ylab
		if (length(ln1) == 1) {
			ln2 <- rep(ln2, nly)
		}
	}
	cells <- ncell(x)


	if (gridded | smooth) {
#			if ((ncell(x) * (nlx + nly)) < .maxmemory()) {
		if ((ncell(x) * (nlx + nly)) < 1000000) {
			maxcell <- ncell(x)
		}
		if (smooth) {
			dots <- list(...)
			if (!is.null(dots$col)) { 
				colramp <- grDevices::colorRampPalette(dots$col)
			}
		}
	}

	x <- as.matrix(spatSample(c(x,y), size=maxcell, method="regular", as.raster=FALSE, warn=FALSE))
	# y <- as.matrix(spatSample(y, size=maxcell, method="regular", as.raster=FALSE))

	y <- x[,c((nlx+1):ncol(x))]
	x <- x[,1:nlx]


	if (warn & (NROW(x) < cells)) {
		warn("plot", 'plot used a sample of ', round(100*NROW(x)/cells, 1), '% of the cells. You can use "maxcell" to increase the sample)')
	}

	if (missing(cex)) {
		if (NROW(x) < 100) {
			cex <- 1
		} else if (NROW(x) < 1000) {
			cex <- 0.5
		} else {
			cex <- 0.2
		}
	}

	if (nlx != nly) {
		# recycling
		d <- cbind(as.vector(x), as.vector(y))
		x <- matrix(d[,1], ncol=nl)
		y <- matrix(d[,2], ncol=nl)
		lab <- vector(length=nl)
		lab[] <- ln1
		ln1 <- lab
		lab[] <- ln2
		ln2 <- lab
	}

	if (nl > 1) {

		old.par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old.par))
		px <- trunc(sqrt(nl))
		py <- ceiling(nl / px)
		graphics::par(mfrow=c(px, py), mar=c(3, 3, 1, 1))

		if (smooth) {
			for (i in 1:nl) {
				graphics::smoothScatter(x[,i], y[,i], main=main[i], xlab=ln1[i], ylab=ln2[i], add=add, asp=asp, colramp=colramp, ...)
			}				
		} else if (gridded) {
			for (i in 1:nl) {
				.plotdens(x[,i], y[,i], nc=ncol, nr=nrow, main=main[i], xlab=ln1[i], ylab=ln2[i], add=add, asp=asp, ...)
			}		
		} else {
			if (add) {
				for (i in 1:nl) {
					graphics::points(x[,i], y[,i], cex=cex, ...)
				}
			} else {
				for (i in 1:nl) {
					plot(x[,i], y[,i], cex=cex, xlab=ln1[i], ylab=ln2[i], main=main[i], asp=asp, ...)
				}
			}
		} 
	} else  {
		if (smooth) {
			graphics::smoothScatter(x, y, main=main[1], xlab=ln1[1], ylab=ln2[1], add=add, asp=asp, colramp=colramp, ...)
		} else if (gridded) {
			.plotdens(x, y, nc=ncol, nr=nrow, main=main[1], xlab=ln1[1], ylab=ln2[1], add=add, asp=asp, ...)
		} else {
			if (add) {
				graphics::points(x, y, cex=cex, ...)
			} else {
				plot(x, y, cex=cex, xlab=ln1[1], ylab=ln2[1], main=main[1], asp=asp, ...)
			}
		}
	}
}


setMethod("plot", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, maxcell=100000, warn=TRUE, nc, nr, maxnl=16, smooth=FALSE, gridded=FALSE, ncol=25, nrow=25, ...) {

		nl <- max(nlyr(x), nlyr(y))
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


		.scatterPlotRaster(x, y, maxcell=maxcell, warn=warn, nc=nc, nr=nr, maxnl=maxnl, gridded=gridded, smooth=smooth, ncol=ncol, nrow=nrow, ...)
	}
)



.plotdens <- function(x, y, nc, nr, xlim=NULL, ylim=NULL, asp=NULL, ...) {
	xy <- stats::na.omit(cbind(x,y))
	if (nrow(xy) == 0) {
		error("plot (density)", "only NA values (in this sample?)")
	}
	r <- apply(xy, 2, range)
	rx <- r[,1]
	if (rx[1] == rx[2]) {
		rx[1] <- rx[1] - 0.5
		rx[2] <- rx[2] + 0.5
	}
	ry <- r[,2]
	if (ry[1] == ry[2]) {
		ry[1] <- ry[1] - 0.5
		ry[2] <- ry[2] + 0.5
	}

	out <- rast(xmin=rx[1], xmax=rx[2], ymin=ry[1], ymax=ry[2], ncol=nc, nrow=nr, crs="+proj=utm +zone=1 +datum=WGS84")
	colnames(xy) <- c("x", "y")
	out <- rasterize(vect(xy), out, fun=function(x, ...) length(x), background=NA)
	if (!is.null(xlim) | !is.null(ylim)) {
		if (is.null(xlim)) xlim <- c(xmin(x), xmax(x))
		if (is.null(ylim)) ylim <- c(ymin(x), ymax(x))
		e <- extent(xlim, ylim)
		out <- extend(crop(out, e), e, value=0)
	}
	plot(out, maxcell=nc*nr, asp=asp, ...)
}


