# comment

setMethod("animate", signature(x="SpatRaster"),
function(x, pause=0.25, main, range=NULL, maxcell=50000, n=1, ...) {
	if (missing(main)) {
		main <- names(x)
	}

#	x <- spatSample(x, size=maxcell, method="regular", as.raster=TRUE, warn=FALSE)
	x <- sampleRaster(x, maxcell, method="regular", replace=FALSE, ext=NULL, warn=FALSE, overview=TRUE)

    if (is.null(range)) {
      range <- range(minmax(x, compute=TRUE))
    } else if (any(is.na(range))) {
      range <- NULL
    }

	nl <- nlyr(x)
	n <- max(1, round(n))
	i <- 1
	reps <- 0
    while (reps < n) {
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

setMethod("animate", signature(x="SpatVector"),
          function(x, pause=0.25, main="", n=1, add=NULL, ...) {
            
            nr <- nrow(x)
            n <- max(1, round(n))
            i <- 1
            reps <- 0
            
            while (reps < n) {
              addd <- if (is.null(add)) i != 1 else add 
              plot(x[i, ], ext=ext(x), add=addd, ...)
              grDevices::dev.flush()
              grDevices::dev.hold()
              Sys.sleep(pause)
              i <- i + 1
              if (i > nr) {
                i <- 1
                reps <- reps+1
              }
            }
          }
)
