# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# Licence GPL v3

setMethod("RGB<-", signature(x="SpatRaster"), 
	function(x, value) {
		if (is.null(value[1]) || is.na(value[1])) {
			x@ptr$removeRGB()
		} else {
			stopifnot(length(value) == 3)
			stopifnot(all(value %in% 1:nlyr(x)))
			
			x@ptr$setRGB(value[1]-1, value[2]-1, value[3]-1)
		}
		messages(x, "RGB<-")
	}
)

setMethod("RGB", signature(x="SpatRaster"), 
	function(x) {
		if (x@ptr$rgb) {
			x@ptr$getRGB() + 1
		} else {
			return(NULL)
		}
	}
)


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
function(x, r=1, g=2, b=3, scale, maxcell=500000, mar=0, stretch=NULL, ext=NULL, smooth=FALSE, colNA="white", alpha, bgalpha, addfun=NULL, zlim=NULL, zlimcol=NULL, axes=FALSE, xlab="", ylab="", asp=NULL, add=FALSE, interpolate, ...) { 

	if (!is.null(mar)) {
		mar <- rep_len(mar, 4)
		if (!any(is.na(mar))) {	
			graphics::par(mar=mar)
		}
	}

	if (missing(scale)) {
		scale <- 255
		if ( all(.hasMinMax(x)) ) {
			rng <- minmax(x)
			scale <- max(max(rng[2]), 255)
		}
	}
	scale <- as.vector(scale)[1]

	if (!is.null(ext)) {
		x <- crop(x, ext)
	}
	x <- spatSample(x[[c(r, g, b)]], maxcell, method="regular", as.raster=TRUE)

	RGB <- values(x)

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
			warn("plotRGB", 'invalid stretch value')
		}
	}


	naind <- as.vector( attr(RGB, "na.action") )
	if (!is.null(naind)) {
		bg <- grDevices::col2rgb(colNA)
		bg <- grDevices::rgb(bg[1], bg[2], bg[3], alpha=bgalpha, max=255)
		z <- rep( bg, times=ncell(x))
		z[-naind] <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, max=scale)
	} else {
		z <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, max=scale)
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


rgb2col <- function(x, r=1, g=2, b=3, filename="", ...) { 
	opt <- spatOptions(filename, ...)
	x@ptr <- x@ptr$rgb2col(r-1, g-1, b-1, opt)
	messages(x, "rgb2col")
}





#### RGB2col

make_cut <- function(x) {
	j <- length(x)
	out <- vector("list", 2*j)
	for (i in 1:j) {
		rgb <- x[[i]]
		if (NROW(rgb) <= 1) {
			out[[i]] <- rgb
			j <- j - 1
			next
		}
		rng <- apply(rgb[,-1], 2, function(i) diff(range(i)))
		if (max(rng) == 0) {
			out[[i]] <- rgb
			j <- j - 1
			next		
		}
		p <- which.max(rng) + 1
		m <- median(rgb[,p])
		out[[i]] <- rgb[rgb[,p] >= m, ,drop=FALSE]
		out[[i+j]] <- rgb[rgb[,p] < m, ,drop=FALSE]
	}
	i <- sapply(out, is.null)
	out <- out[!i]
	i <- sapply(out, nrow) > 0
	out[i]
}

median_cut <- function(v) {
	v <- list(v)
	n <- 0
	while ((length(v) < 129) & (length(v) > n)) {
		n <- length(v)
		v <- make_cut(v)
	}
	s <- sapply(v, function(i) max(apply(i[,-1,drop=FALSE], 2, function(j) diff(range(j)))))
	n <- 256 - length(v)
	ss <- rev(sort(s))
	ss <- max(2, min(ss[1:n]))
	i <- which(s > ss)
	if (length(i) > 0) {
		vv <- make_cut(v[i])
		v <- c(v[-i], vv)
	}
	v <- lapply(1:length(v), function(i) cbind(group=i, v[[i]]))
	do.call(rbind, v)
}


setMethod("RGB2col", signature(x="SpatRaster"), 
	function(x, value, stretch=NULL, grays=FALSE, filename="", overwrite=FALSE, ...) {
		idx <- RGB(x)
		if (is.null(idx)) {
			if (missing(value)) {
				error("rgb2coltab", "x does not have an RGB attribute and the value argument is missing")
			} else {
				idx <- value
			}
		}
		stopifnot(length(idx) == 3)
		if ((min(idx) < 1) | (max(idx) > nlyr(x))) {
			error("rgb2coltab", "invalid value (RGB indices)")	
		}
		x <- x[[idx]]
		
		if (!is.null(stretch)) {
			stretch = tolower(stretch)
			if (stretch == "lin") {
				values(x[[1]]) <- .linStretch(values(x[[1]]))
				values(x[[2]]) <- .linStretch(values(x[[2]]))
				values(x[[3]]) <- .linStretch(values(x[[3]]))
				scale <- 255
			} else if (stretch == "hist") {
				values(x[[1]]) <- .eqStretch(values(x[[1]]))
				values(x[[2]]) <- .eqStretch(values(x[[2]]))
				values(x[[3]]) <- .eqStretch(values(x[[3]]))
				scale <- 255
			} else if (stretch != "") {
				warn("plotRGB", 'invalid stretch value')
			}
		}
		
		
		if (grays) {
			opt <- spatOptions(filename, overwrite, ...)
			x@ptr <- x@ptr$rgb2col(0, 1, 2, opt)
			return(x)
		}
		
		v <- cbind(id=1:ncell(x), values(x))
		v <- median_cut(stats::na.omit(v))
		
		a <- aggregate(v[,3:5], list(v[,1]), median)
		a$cols <- grDevices::rgb(a[,2], a[,3], a[,4], max=255)
		m <- merge(v[,1:2], a[, c(1,5)], by=1)
		r <- rast(x, 1)
		r[m$id] <- m$group - 1
		coltab(r) <- a$cols
		if (filename != "") {
			r <- writeRaster(r, filename, overwrite, ...)
		}
		r
	}
)


