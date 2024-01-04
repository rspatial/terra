# Author: Robert J. Hijmans
# Date :  June 2008
# Version 1.0
# License GPL v3


setMethod("hist", signature(x="SpatRaster"),
	function(x, layer, maxcell=1000000, plot=TRUE, maxnl=16, main, ...) {

		if (missing(layer)) {
			y <- 1:nlyr(x)
		} else if (is.character(layer)) {
			y <- match(layer, names(x))
			maxnl <- Inf
		} else {
			y <- layer
			maxnl <- Inf
		}

		y <- as.integer(round(y))
		y <- stats::na.omit(y)
		y <- y[ y >= 1 & y <= nlyr(x) ]
		nl <- length(y)

		if (nl == 0) {
			error("hist", "no valid layers selected")
		}

		if (missing(main)) {
			main=names(x)
		}

		if (nl > 1)	{
			res <- list()
			if (nl > maxnl) {
				warn(paste("hist", "only the first", maxnl, "layers are used (see argument maxnl)"))
				nl <- maxnl
				y <- y[1:maxnl]
			}

			nc <- ceiling(sqrt(nl))
			nr <- ceiling(nl / nc)
			mfrow <- graphics::par("mfrow")
			spots <- mfrow[1] * mfrow[2]
			if (spots < nl) {
				old.par <- graphics::par(no.readonly =TRUE)
				on.exit(graphics::par(old.par))
				graphics::par(mfrow=c(nr, nc))
			}
			for (i in 1:length(y)) {
				res[[i]] = .hist1(x[[ y[i] ]], maxcell=maxcell, main=main[y[i]], plot=plot, ...)
			}

		} else if (nl==1) {
			if (nlyr(x) > 1) {
				x <- x[[y]]
				main <- main[y]
			}
			res <- .hist1(x, maxcell=maxcell, main=main, plot=plot, ...)
		}
		if (plot) {
			return(invisible(res))
		} else {
			return(res)
		}
	}
)



.hist1 <- function(x, maxcell, main, plot, ...){

	stopifnot(hasValues(x))

	if ( ncell(x) <= maxcell ) {
		v <- stats::na.omit(values(x))
	} else {
		# TO DO: make a function that does this by block and combines  all data into a single histogram
		v <- spatSample(x, maxcell, method="regular", as.df=FALSE, as.raster=FALSE, warn=FALSE)
		msg <- paste("a sample of ", round(100 * length(v) / ncell(x)), "% of the cells was used", sep="")
		if (any(is.na(v))) {
			v <- stats::na.omit(v)
			msg <- paste(msg, " (of which ", 100 - round(100 * length(v) / maxcell ), "% was NA)", sep="")
		}
		warn("hist", msg)
	}

	if (nrow(v) == 0) {
		error("hist", "all (sampled) values are NA")
	}

#	if (.shortDataType(x) == 'LOG') {
#		v <- v * 1
#	}

	if (plot) {
		out <- hist(v, main=main, plot=plot, ...)
	} else {
		out <- hist(v, plot=plot, ...)
	}
	out$xname <- names(x)
	out
}



