

setMethod("animate", signature(x="SpatRaster"), 
function(x, pause=0.25, main, range, maxcell=50000, n=10, ...) {
	nl <- nlyr(x)
	if (missing(main)) {
		main <- names(x)
	}

	x <- spatSample(x, size=maxcell, method="regular", as.raster=TRUE)
	
	if (missing(range)) {
		mnmx <- minmax(x)
		range <- c(min(mnmx[1,]), max(mnmx[2,]))
	}
	
	i <- 1
	reps <- 0
    while (reps < n) {
        plot(x[[i]], main = main[i], range=range, maxcell=Inf, ...)
        grDevices::dev.flush()
        Sys.sleep(pause)
        i <- i + 1
        if (i > nl) {
            i <- 1
			reps <- reps+1
		}
    }
}
)
