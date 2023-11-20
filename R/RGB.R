# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# License GPL v3


setMethod ("has.RGB" , "SpatRaster",
	function(x, strict=TRUE) {
		if (strict) {
			x@cpp$rgbtype == "rgb"
		} else {
			x@cpp$rgbtype != ""
		}
	}
)


setMethod("set.RGB", signature(x="SpatRaster"),
	function(x, value=1:3, type="rgb") {
		if (is.null(value[1]) || is.na(value[1]) || any(value < 1)) {
			x@cpp$removeRGB()
		} else {
			stopifnot(all(value %in% 1:nlyr(x)))
			if (length(value) == 3) {
				x@cpp$setRGB(value[1]-1, value[2]-1, value[3]-1, -99, type)
			} else if (length(value) == 4) {
				x@cpp$setRGB(value[1]-1, value[2]-1, value[3]-1, value[4]-1, type)
			} else {
				error("set.RGB", "value must have length 3 or 4")
			}
		}
		x <- messages(x, "set.RGB")
		invisible(TRUE)
	}
)

setMethod("RGB<-", signature(x="SpatRaster"),
	function(x, ..., type="rgb", value) {
		x@cpp <- x@cpp$deepcopy()
		set.RGB(x, value, type)
		x
	}
)

setMethod("RGB", signature(x="SpatRaster"),
	function(x) {
		if (x@cpp$rgb) {
			x@cpp$getRGB() + 1
		} else {
			return(NULL)
		}
	}
)

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


rgb2col <- function(x, value, stretch=NULL, grays=FALSE, NAzero=FALSE, filename="", overwrite=FALSE, ...) {
	idx <- RGB(x)
	if (is.null(idx)) {
		if (missing(value)) {
			error("colorize", "x does not have an RGB attribute and the value argument is missing")
		} else {
			idx <- value
		}
	}
	n <- length(idx)
	stopifnot((n == 3) | (n == 4))
	if ((min(idx) < 1) | (max(idx) > nlyr(x))) {
		error("colorize", "invalid value (RGB indices)")
	}
	x <- x[[idx]]

	if (!is.null(stretch)) {
		if (stretch == "lin") {
			x <- stretch(x, minq=0.02, maxq=0.98)
		} else {
			x <- stretch(x, histeq=TRUE, scale=255)
		}
	}

	if (n == 4) x[[4]] <- x[[4]] * 255

	if (NAzero) {
		x <- classify(x, cbind(NA, 0))
	}


	if (grays) {
		opt <- spatOptions(filename, overwrite, ...)
		x@cpp <- x@cpp$rgb2col(0, 1, 2, opt)
		return(messages(x, "colorize"))
	}

	v <- cbind(id=1:ncell(x), values(x))
	v <- median_cut(stats::na.omit(v))

	a <- aggregate(v[,-c(1:2)], list(v[,1]), median)
#	if (n==3) {
#		a$cols <- grDevices::rgb(a[,2], a[,3], a[,4], maxColorValue=255)
#	} else {
#		a$cols <- grDevices::rgb(a[,2], a[,3], a[,4], a[,5], maxColorValue=255)
#	}
	m <- merge(v[,1:2], a, by=1)
	m[,1] <- m[,1] - 1
	r <- rast(x, 1)
	r[m$id] <- m$group
	coltab(r) <- m[,-2]
	if (filename != "") {
		r <- writeRaster(r, filename, overwrite, ...)
	}
	r
}


terra_col2rgb <- function(x, alpha=FALSE, filename="", overwrite=FALSE, ...) {
	if (nlyr(x) > 1) {
		x <- x[[1]]
		warn("colorize", "only the first layer of 'x' is considered")
	}
	ct <- coltab(x)[[1]]
	if (is.null(ct)) {
		error("error", "x has no color table")
	}
	ct <- as.matrix(ct)
	nms <- c("red", "green", "blue", "alpha")
	rgbidx <- 1:4
	if (!alpha) {
		ct <- ct[,1:4]
		nms <- nms[1:3]
		rgbidx <- rgbidx[1:3]
	}

	wopt=list(...)
	if (is.null(wopt$names)) {
		wopt$names <- nms
	}
	out <- subst(x, from=ct[,1], to=ct[,-1], raw=TRUE, filename=filename, overwrite=overwrite, wopt=wopt)
	set.RGB(out, rgbidx)
	out
}



setMethod("colorize", signature(x="SpatRaster"),
	function(x, to="hsv", alpha=FALSE, stretch=NULL, grays=FALSE, NAzero=FALSE, filename="", overwrite=FALSE, ...) {
		to <- tolower(to)
		if (to %in% c("hsi", "hsl", "hsv")) {
			opt <- spatOptions(filename, overwrite, ...)
			x@cpp <- x@cpp$rgb2hsx(to, opt)
		} else if (to == "rgb") {
			if (nlyr(x) == 1) {
				return(terra_col2rgb(x, alpha=alpha, filename=filename, overwrite=overwrite, ...))
			} else {
				opt <- spatOptions(filename, overwrite, ...)
				x@cpp <- x@cpp$hsx2rgb(opt)
			}
		} else if (to == "hsl") {
			opt <- spatOptions(filename, overwrite, ...)
			x@cpp <- x@cpp$hsx2rgb(to, opt)
		} else if (to == "col") {
			return(rgb2col(x, stretch=stretch, grays=grays, NAzero=NAzero, filename=filename, overwrite=overwrite, ...))
		}
		messages(x)
	}
)

