# Author: Robert J. Hijmans
# Date : June 2013
# Version 1.0
# License GPL v3


.circular.weight <- function(rs, d, fillNA=FALSE) {
	nx <- 1 + 2 * floor(d/rs[1])
	ny <- 1 + 2 * floor(d/rs[2])
	m <- matrix(ncol=nx, nrow=ny)
	m[ceiling(ny/2), ceiling(nx/2)] <- 1
	if ((nx != 1) || (ny != 1)) {
		x <- rast(m, crs="+proj=utm +zone=1 +datum=WGS84")
		ext(x) <- c(xmin=0, xmax=nx*rs[1], ymin=0, ymax=ny*rs[2])
		d <- as.matrix(distance(x), wide=TRUE) <= d
		m <- d / sum(d)
	}
	if (fillNA) {
		m[m <= 0] <- NA
	}
	m
}



.Gauss.weight <- function(rs, sigma) {
	if (length(sigma) == 1) {
		d <- 3 * sigma
	} else {
		d <- sigma[2]
		sigma <- sigma[1]
	}
	nx <- 1 + 2 * floor(d/rs[1])
	ny <- 1 + 2 * floor(d/rs[2])
	m <- matrix(ncol=nx, nrow=ny)
	xr <- (nx * rs[1]) / 2
	yr <- (ny * rs[2]) / 2
	r <- rast(m, crs="+proj=utm +zone=1 +datum=WGS84")
	ext(r) <- c(xmin=-xr[1], xmax=xr[1], ymin=-yr[1], ymax=yr[1])
	p <- xyFromCell(r, 1:ncell(r))^2
# according to http://en.wikipedia.org/wiki/Gaussian_filter
	m <- 1/(2*pi*sigma^2) * exp(-(p[,1]+p[,2])/(2*sigma^2))
	m <- matrix(m, ncol=nx, nrow=ny, byrow=TRUE)
# sum of weights should add up to 1
	m / sum(m)
}


.rectangle.weight <- function(rs, d) {
	d <- rep(d, length.out=2)
	nx <- 1 + 2 * floor(d[1]/rs[1])
	ny <- 1 + 2 * floor(d[2]/rs[2])
	m <- matrix(1, ncol=nx, nrow=ny)
	m / sum(m)
}



focalMat <- function(x, d, type=c('circle', 'Gauss', 'rectangle'), fillNA=FALSE) {
	type <- match.arg(type)
	x <- res(x)
	if (type == 'circle') {
		.circular.weight(x, d[1], fillNA)
	} else if (type == 'Gauss') {
		if (!length(d) %in% 1:2) {
			error("focalMath", "if type=Gauss, d should be a vector of length 1 or 2")
		}
		.Gauss.weight(x, d)
	} else {
		.rectangle.weight(x, d)
	}
}


