# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# Licence GPL v3


.linStretch <- function (x) {
    v <- stats::quantile(x, c(0.02, 0.98), na.rm = TRUE)
    temp <- (255 * (x - v[1]))/(v[2] - v[1])
    temp[temp < 0] <- 0
    temp[temp > 255] <- 255
    return(temp)
}

# Histogram equalization stretch
.eqStretch <- function(x){
	ecdfun <- stats::ecdf(x)
	ecdfun(x)*255
}



setMethod("plotRGB", signature(x="SpatRaster"), 
function(x, r=1, g=2, b=3, scale, maxcell=500000, stretch=NULL, ext=NULL, interpolate=FALSE, colNA='white', alpha, bgalpha, addfun=NULL, zlim=NULL, zlimcol=NULL, axes=FALSE, xlab='', ylab='', asp=NULL, add=FALSE, margins=FALSE, ...) { 

	if (missing(scale)) {
		scale <- 255
		if ( all(.hasMinMax(x)) ) {
			rng <- minmax(x)
			scale <- max(max(rng[2]), 255)
		}
	}
	scale <- as.vector(scale)[1]
	
	r <- spatSample(x[[r]], maxcell, ext=ext, method="regular", as.raster=TRUE)
	g <- spatSample(x[[g]], maxcell, ext=ext, method="regular", as.raster=TRUE)
	b <- spatSample(x[[b]], maxcell, ext=ext, method="regular", as.raster=TRUE)

	RGB <- cbind(values(r), values(g), values(b))
	
	if (!is.null(zlim)) {
		if (length(zlim) == 2) {
			zlim <- sort(zlim)
			if (is.null(zlimcol)) {
				RGB[ RGB<zlim[1] ] <- zlim[1]
				RGB[ RGB>zlim[2] ] <- zlim[2]
			} else { #if (is.na(zlimcol)) {
				RGB[RGB<zlim[1] | RGB>zlim[2]] <- NA
			} 
		} else if (NROW(zlim) == 3 & NCOL(zlim) == 2) {
			for (i in 1:3) {
				zmin <- min(zlim[i,])		
				zmax <- max(zlim[i,])
				if (is.null(zlimcol)) {
					RGB[RGB[,i] < zmin, i] <- zmin
					RGB[RGB[,i] > zmax, i] <- zmax
				} else { #if (is.na(zlimcol)) {
					RGB[RGB < zmin | RGB > zmax, i] <- NA
				}
			}
		} else {
			stop('zlim should be a vector of two numbers or a 3x2 matrix (one row for each color)')
		}
	}
	
	RGB <- stats::na.omit(RGB)
	
	if (!is.null(stretch)) {
		stretch = tolower(stretch)
		if (stretch == 'lin') {
			RGB[,1] <- .linStretch(RGB[,1])
			RGB[,2] <- .linStretch(RGB[,2])
			RGB[,3] <- .linStretch(RGB[,3])
			scale <- 255
		} else if (stretch == 'hist') {
			RGB[,1] <- .eqStretch(RGB[,1])
			RGB[,2] <- .eqStretch(RGB[,2])
			RGB[,3] <- .eqStretch(RGB[,3])
			scale <- 255
		} else if (stretch != '') {
			warning('invalid stretch value')
		}
	}

	
	naind <- as.vector( attr(RGB, "na.action") )
	if (!is.null(naind)) {
		bg <- grDevices::col2rgb(colNA)
		bg <- grDevices::rgb(bg[1], bg[2], bg[3], alpha=bgalpha, max=255)
		z <- rep( bg, times=ncell(r))
		z[-naind] <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, max=scale)
	} else {
		z <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, max=scale)
	}
	
	z <- matrix(z, nrow=nrow(r), ncol=ncol(r), byrow=T)

	requireNamespace("grDevices")
	bb <- as.vector(matrix(as.vector(ext(x)), ncol=2))

	bb <- as.vector(ext(r))
	
	if (!add) {
		if ((!axes) & (!margins)) {
			graphics::par(plt=c(0,1,0,1))
		}

		if (is.null(asp)) {
			if (couldBeLonLat(x)) {
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
			if (xres(r) %% 1 == 0) xticks = round(xticks)
			if (yres(r) %% 1 == 0) yticks = round(yticks)
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
)

