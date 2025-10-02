# comment

setMethod("animate", signature(x="SpatRaster"),
function(x, pause=0.25, main, range, maxcell=50000, n=1, ...) {
	if (missing(main)) {
		main <- names(x)
	}

#	x <- spatSample(x, size=maxcell, method="regular", as.raster=TRUE, warn=FALSE)
	x <- sampleRaster(x, maxcell, method="regular", replace=FALSE, ext=NULL, warn=FALSE, overview=TRUE)

    if (missing(range)) {
      mnmx <- minmax(x)
      range <- c(min(mnmx[1,]), max(mnmx[2,]))
    }
    if (length(range) == 1 && is.na(range)) {
      range <- minmax(x)
    }

	nl <- nlyr(x)
	n <- max(1, round(n))
	i <- 1
	reps <- 0
    while (reps < n) {
        if (is.matrix(range))
          plot(x[[i]], main = main[i], range=range[, i], maxcell=Inf, ...)
        else 
          plot(x[[i]], main = main[i], range=range, maxcell=Inf, ...)
        grDevices::dev.flush()
	    grDevices::dev.hold()
        Sys.sleep(pause)
        i <- i + 1
        if (i > nl) {
            i <- 1
			reps <- reps+1
		}
    }
}
)
