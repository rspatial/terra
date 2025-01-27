# Author: Robert J. Hijmans
# Date : January 2025
# Version 1.0
# License GPL v3


otsu <- function(values, counts) {
#based on https://stackoverflow.com/questions/51116495/auto-thresholding-on-r-raster-object by Eric
	n <- length(values)
	w1 <- cumsum(counts)
	w2 <- w1[n] + counts - w1
	cv <- counts * values
	m1 <- cumsum(cv)
	m2 <- m1[n] + cv - m1
	varian <- w1 * w2 * (m2/w2 - m1/w1)^2
	mxi <- which(varian == max(varian, na.rm = TRUE))
	(values[mxi[1]] + values[mxi[length(mxi)]])/2
}


setMethod("twoClasses", signature(x="SpatRaster"),
	function(x, method="otsu", maxcell=1000000, combine=FALSE, as.raster=TRUE, filename="", ...) {

		method <- match.arg(tolower(method), c("mean", "median", "otsu"))
		nl <- nlyr(x)
		if (combine) {
			if (method == "otsu") {
				rng <- minmax(x, compute=TRUE)
				breaks <- seq(rng[1], rng[2], length.out=257)
				h <- classify(x, breaks, include.lowest=TRUE, right=FALSE)
				f <- freq(h)
				f <- aggregate(f[,3,drop=FALSE], f[,2,drop=FALSE], sum)
				breaks <- breaks + ((breaks[2] - breaks[1]) / 2)
				f$value <- breaks[f$value+1]
				th <- otsu(f$value, f$count)
			} else if (method == "mean") {
				r <- spatSample(x, maxcell, "regular", as.df=FALSE)
				th <- mean(r, na.rm=TRUE)
			} else if (method == "median") {
				r <- spatSample(x, maxcell, "regular", as.df=FALSE)
				th <- median(r, na.rm=TRUE)
			}
		} else {
			if (method == "otsu") {
				th <- rep(NA, nl)
				rng <- minmax(x, compute=TRUE)
				for (i in 1:nl) {
					breaks <- seq(rng[1,i], rng[2, i], length.out=257)
					h <- classify(x, breaks, include.lowest=TRUE, right=FALSE)
					f <- freq(h)
					breaks <- breaks + ((breaks[2] - breaks[1]) / 2)
					f[,2] <- breaks[f[,2]+1]
					j <- f[,1] == i
					th[i] <- otsu(f[j,2], f[j,3])
				} 
			} else if (method == "mean") {
				th <- global(x, "mean", na.rm=TRUE)[,1]
			} else if (method == "median") {
				r <- spatSample(x, maxcell, "regular", as.df=FALSE)
				th <- apply(r, 2, median, na.rm=TRUE)
			}
		}
		
		if (!as.raster) {
			if (combine) {
				th <- as.vector(th)
			} else {
				names(th) <- names(x)
			}
			return(th)
		} 
		opt <- spatOptions(filename, ...)
		x@pntr <- x@pntr$arith_numb(th, ">", FALSE, FALSE, opt)
		messages(x, "twoClasses")
	}
)
