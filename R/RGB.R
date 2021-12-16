# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# License GPL v3




setMethod("RGB<-", signature(x="SpatRaster"), 
	function(x, value) {
		if (is.null(value[1]) || is.na(value[1])) {
			x@ptr$removeRGB()
		} else {
			stopifnot(all(value %in% 1:nlyr(x)))
			if (length(value) == 3) {
				x@ptr$setRGB(value[1]-1, value[2]-1, value[3]-1, -99)
			} else if (length(value) == 4) {
				x@ptr$setRGB(value[1]-1, value[2]-1, value[3]-1, value[4]-1)
			} else {
				error("RGB<-", "value must have length 3 or 4")
			}
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
		n <- length(idx)
		stopifnot((n == 3) | (n == 4))
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

		if (n == 4) x[[4]] <- x[[4]] * 255

		if (grays) {
			opt <- spatOptions(filename, overwrite, ...)
			x@ptr <- x@ptr$rgb2col(0, 1, 2, opt)
			return(messages(x, "RGB2col"))
		}

		v <- cbind(id=1:ncell(x), values(x))
		v <- median_cut(stats::na.omit(v))

		a <- aggregate(v[,-c(1:2)], list(v[,1]), median)
		if (n==3) {
			a$cols <- grDevices::rgb(a[,2], a[,3], a[,4], maxColorValue=255)
		} else {
			a$cols <- grDevices::rgb(a[,2], a[,3], a[,4], a[,5], maxColorValue=255)
		}
		m <- merge(v[,1:2], a[, c(1,n+2)], by=1)
		r <- rast(x, 1)
		r[m$id] <- m$group - 1
		coltab(r) <- a$cols
		if (filename != "") {
			r <- writeRaster(r, filename, overwrite, ...)
		}
		r
	}
)


setMethod("col2RGB", signature(x="SpatRaster"), 
	function(x, alpha=FALSE, filename="", overwrite=FALSE, ...) {
		if (nlyr(x) > 1) {
			x <- x[[1]]
			warn("col2RGB", "only the first layer of 'x' is considered")
		}
		ct <- coltab(r)[[1]]
		if (is.null(ct)) {
			error("error", "x has no color table")
		}
		ct <- as.matrix(ct)
		if (!alpha) {
			ct <- ct[,1:3]
		}
		r <- app(x, function(i) { ct[i+1, ,drop=FALSE] }, filename="", overwrite=FALSE, wopt=list(...))
		RGB(r) <- 1:3
		r
	}
)


setMethod("RGB2HS", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, ...) {
	}
)

setMethod("HS2RGB", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, ...) {
	}
)

