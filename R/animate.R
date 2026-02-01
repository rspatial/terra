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

	n <- max(1, round(n))
    for (reps in 1:n) {
		for (i in 1:nlyr(x)) {
			plot(x[[i]], main = main[i], range=range, maxcell=Inf, ...)
			grDevices::dev.flush()
			grDevices::dev.hold()
			Sys.sleep(pause)
		}
	}
	g <- grDevices::dev.flush()
}
)

setMethod("animate", signature(x="SpatVector"),
function(x, pause=0.25, main, n=1, vars=NULL, range=NULL, add=NULL, ...) {

	n <- max(1, round(n))
		
	if (!is.null(vars)) {

		if (is.null(add)) add <- TRUE
		if (is.numeric(vars)) {
			vars <- names(x)[vars]
		} else {
			stopifnot(vars %in% names(x))
		}

		if (missing(main)) {
			main <- vars
		} else {
			main <- rep_len(main, length(vars))
		}

		if (is.null(range)) {
			v <- values(x[,vars])
			d <- sapply(v, is.numeric)
			if (all(d)) {
				range <- range(v)
			} 
		} else if (any(is.na(range))) {
			range <- NULL
		}

		for (reps in 1:n) {
			for (i in 1:length(vars)) {
				plot(x, vars[i], main=main[i], add=add)
				grDevices::dev.flush()
				grDevices::dev.hold()
				Sys.sleep(pause)
			}
		}
		
	} else {

		if (missing(main)) {
			main <- ""
		}

		for (reps in 1:n) {
			for (i in 1:nrow(x)) {
				addd <- if (is.null(add)) i != 1 else add 
				plot(x[i, ], ext=ext(x), add=addd, ...)
				grDevices::dev.flush()
				grDevices::dev.hold()
				Sys.sleep(pause)
			}
		}
	}
	g <- grDevices::dev.flush()
}
)
