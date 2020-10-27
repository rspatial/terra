# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# Licence GPL v3


.plotCT <- function(x, maxcell=500000, stretch=NULL, ext=NULL, interpolate=FALSE, colNA='white', alpha, bgalpha, addfun=NULL, zlim=NULL, zlimcol=NULL, axes=FALSE, xlab='', ylab='', asp=NULL, add=FALSE, margins=FALSE, ...) { 


	d <- x@ptr$getColors()
	d <- terra:::.getSpatDF(d[[1]])
	ct <- grDevices::rgb(d[,1], d[,2], d[,3], d[,4], maxColorValue=255)

	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	x <- spatSample(x[[1]], maxcell, method="regular", as.raster=TRUE)

	z <- as.integer(values(x))
	z[z<0 | z>255] <- NA
	
	z <- ct[z]
	z <- matrix(z, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)

	requireNamespace("grDevices")
	bb <- as.vector(matrix(as.vector(ext(x)), ncol=2))

	bb <- as.vector(ext(x))
	
	if (!add) {
		#if ((!axes) & (!margins)) {
		#	old.par <- graphics::par(no.readonly =TRUE)
		#	on.exit(graphics::par(old.par))   
		#	graphics::par(plt=c(0,1,0,1))
		#}

		if (is.null(asp)) {
			if (isLonLat(x, perhaps=TRUE, warn=FALSE)) {
			    ym <- mean(bb[3:4])
				asp <- 1/cos((ym * pi)/180)
				#asp <- min(5, 1/cos((ym * pi)/180))
			} else {
				asp <- 1
			}
		}
		
		xlim=c(bb[1], bb[2])
		ylim=c(bb[3], bb[4])
		
		plot(NA, NA, xlim=xlim, ylim=ylim, type = "n", xaxs='i', yaxs='i', xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, ...)
		if (axes) {
			xticks <- graphics::axTicks(1, c(xlim[1], xlim[2], 4))
			yticks <- graphics::axTicks(2, c(ylim[1], ylim[2], 4))
			
			if (xres(x) %% 1 == 0) xticks = round(xticks)
			if (yres(x) %% 1 == 0) yticks = round(yticks)
			graphics::axis(1, at=xticks)
			graphics::axis(2, at=yticks, las = 1)
			#graphics::axis(3, at=xticks, labels=FALSE, lwd.ticks=0)
			#graphics::axis(4, at=yticks, labels=FALSE, lwd.ticks=0)
		}
	}
	graphics::rasterImage(z, bb[1], bb[3], bb[2], bb[4], interpolate=interpolate, ...)
	
	if (!is.null(addfun)) {
		if (is.function(addfun)) {
			addfun()
		}
	}
}


