


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

rgb2coltab <- function(x, value, grays=FALSE, filename="", overwrite=FALSE, ...) {
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
	
	if (grays) {
		opt <- terra:::spatOptions(filename, overwrite, ...)
		x@ptr <- x@ptr$rgb2col(0, 1, 2, opt)
		return(x)
	}
	
	v <- cbind(id=1:ncell(x), values(x))
	v <- median_cut(na.omit(v))
	
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

