# Author: Robert J. Hijmans
# Date : June 2013
# Version 1.0
# Licence GPL v3


.circular.weight <- function(rs, d) {
	nx <- 1 + 2 * floor(d/rs[1])
	ny <- 1 + 2 * floor(d/rs[2])
	m <- matrix(ncol=nx, nrow=ny)
	m[ceiling(ny/2), ceiling(nx/2)] <- 1
	if (nx == 1 & ny == 1) {
		return(m)
	} else {
		x <- raster(m, xmn=0, xmx=nx*rs[1], ymn=0, ymx=ny*rs[2], crs=CRS("+proj=utm +zone=1 +datum=WGS84"))
		d <- as.matrix(distance(x)) <= d
		d / sum(d)
	}
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
	r <- raster(m, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs=CRS("+proj=utm +zone=1 +datum=WGS84"))
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



focalMat <- function(x, d, type=c('circle', 'Gauss', 'rectangle')) {
	type <- match.arg(type)
	x <- res(x)
	if (type == 'circle') {
		.circular.weight(x, d[1])
	} else if (type == 'Gauss') {
		if (!length(d) %in% 1:2) {
			stop("If type=Gauss, d should be a vector of length 1 or 2")
		}
		.Gauss.weight(x, d)
	} else {
		.rectangle.weight(x, d)
	}
}


