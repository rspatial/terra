# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# License GPL v3

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

rgbstretch <- function(RGB, stretch, caller="") {
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
		warn(caller, "invalid stretch value")
	}
}



setMethod("plotRGB", signature(x="SpatRaster"), 
function(x, r=1, g=2, b=3, a=NULL, scale, maxcell=500000, mar=0, stretch=NULL, ext=NULL, smooth=FALSE, colNA="white", alpha, bgalpha, addfun=NULL, zlim=NULL, zlimcol=NULL, axes=FALSE, xlab="", ylab="", asp=NULL, add=FALSE, interpolate, ...) { 

	x <- x[[c(r, g, b, a)]]

	if (!is.null(mar)) {
		mar <- rep_len(mar, 4)
		if (!any(is.na(mar))) {
			graphics::par(mar=mar)
		}
	}

	if (missing(scale)) {
		scale <- 255
		if ( all(hasMinMax(x)) ) {
			rng <- minmax(x)[, 1:3]
			scale <- max(max(rng[2]), 255)
		}
	}
	scale <- as.vector(scale)[1]

	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	x <- spatSample(x, maxcell, method="regular", as.raster=TRUE)

	RGB <- values(x)
	RGB <- stats::na.omit(RGB)
	naind <- as.vector( attr(RGB, "na.action") )

	if (!is.null(a)) {
		alpha <- RGB[,4] * 255
		RGB <- RGB[,-4]
	}

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
			error('zlim should be a vector of two numbers or a 3x2 matrix (one row for each color)')
		}
	}


	if (!is.null(stretch)) {
		RGB <- rgbstretch(RGB, stretch, "plotRGB")
	}


	if (!is.null(naind)) {
		bg <- grDevices::col2rgb(colNA)
		bg <- grDevices::rgb(bg[1], bg[2], bg[3], alpha=bgalpha, maxColorValue=255)
		z <- rep( bg, times=ncell(x))
		z[-naind] <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, maxColorValue=scale)
	} else {
		z <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, maxColorValue=scale)
	}

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
			if (is.lonlat(x, perhaps=TRUE, warn=FALSE)) {
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
	if (!missing(interpolate)) { # for backwards compatibility
		if (is.logical(interpolate)) {
			smooth <- interpolate
		}
	}
	graphics::rasterImage(z, bb[1], bb[3], bb[2], bb[4], interpolate=smooth, ...)

	if (!is.null(addfun)) {
		if (is.function(addfun)) {
			addfun()
		}
	}
}
)

