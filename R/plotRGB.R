# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# License GPL v3


.as.raster.rgb <- function(out, x) {

	if (is.null(out$rgb$scale)) {
		scale <- 255
		if ( all(hasMinMax(x)) ) {
			rng <- minmax(x)[, 1:3]
			scale <- max(max(rng[2]), 255)
		}
	} else {
		scale <- out$rgb$scale
	}
	
	if (!is.null(out$rgb$zlim)) {
		if (length(out$rgb$zlim) == 2) {
			out$rgb$zlim <- sort(out$rgb$zlim)
			if (isTRUE(out$rgb$zcol)) {
				x <- clamp(x, out$rgb$zlim[1], out$rgb$zlim[2], values=TRUE)
			} else { #if (is.na(zlimcol)) {
				x <- clamp(x, out$rgb$zlim[1], out$rgb$zlim[2], values=FALSE)
			}
		} else if (NROW(out$rgb$zlim) == 3 & NCOL(out$rgb$zlim) == 2) {
			for (i in 1:3) {
				zmin <- min(out$rgb$zlim[i,])
				zmax <- max(out$rgb$zlim[i,])
				if (isTRUE(out$rgb$zcol)) {
					x[[i]] <- clamp(x[[i]], zmin, zmax, values=TRUE)
				} else { #if (is.na(zlimcol)) {
					x[[i]] <- clamp(x[[i]], zmin, zmax, values=FALSE)
				}
			}
		} else {
			error('zlim should be a vector of two numbers or a 3x2 matrix (one row for each color)')
		}
	}

	if (!is.null(out$rgb$stretch)) {
		if (out$rgb$stretch == "lin") {
			if ((!is.null(out$rgb$zlim)) && (length(out$rgb$zlim) == 2)) {
				x <- stretch(x, smin=out$rgb$zlim[1], smax=out$rgb$zlim[2])
			} else {
				x <- stretch(x, minq=0.02, maxq=0.98)
			}
		} else {
			x <- stretch(x, histeq=TRUE, scale=255)
		}
		scale <- 255
	}

	RGB <- values(x)
	RGB <- stats::na.omit(RGB)
	naind <- as.vector( attr(RGB, "na.action") )

	if (ncol(RGB) == 4) {
		alpha <- RGB[,4]
		RGB <- RGB[,-4]
	} else {
		alpha <- out$alpha
	}

	if (!is.null(naind)) {
		bg <- grDevices::col2rgb(out$rgb$colNA)
		if (is.null(out$rgb$bgalpha)) out$rgb$bgalpha <- 255
		bg <- grDevices::rgb(bg[1], bg[2], bg[3], alpha=out$rgb$bgalpha, maxColorValue=255)
		z <- rep( bg, times=ncell(x))
		z[-naind] <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, maxColorValue=scale)
	} else {
		z <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, maxColorValue=scale)
	}
	
	out$r <- matrix(z, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
	out
}

# ..linStretch <- function (x) {
    # v <- stats::quantile(x, c(0.02, 0.98), na.rm = TRUE)
    # temp <- (255 * (x - v[1]))/(v[2] - v[1])
    # temp[temp < 0] <- 0
    # temp[temp > 255] <- 255
    # return(temp)
# }

# # Histogram equalization stretch
# ..eqStretch <- function(x){
	# ecdfun <- stats::ecdf(x)
	# ecdfun(x)*255
# }

# ..rgbstretch <- function(RGB, stretch, caller="") {
	# stretch = tolower(stretch)
	# if (stretch == 'lin') {
		# RGB[,1] <- ..linStretch(RGB[,1])
		# RGB[,2] <- ..linStretch(RGB[,2])
		# RGB[,3] <- ..linStretch(RGB[,3])
	# } else if (stretch == 'hist') {
		# RGB[,1] <- ..eqStretch(RGB[,1])
		# RGB[,2] <- ..eqStretch(RGB[,2])
		# RGB[,3] <- ..eqStretch(RGB[,3])
	# } else if (stretch != '') {
		# warn(caller, "invalid stretch value")
	# }
	# RGB
# }


