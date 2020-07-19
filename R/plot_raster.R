

.plot.legend <- function(x) {
	if (is.null(x$leg$ext)) {
		x <- .get.leg.extent(x)
	} else {
		x <- .get.leg.coords(x)	
	}
	if (x$leg$type == "continuous") {
		x <- .plot.cont.legend(x)
	}
	x
}


.plotit <- function(x, leg.ext=NULL, leg.levels=NULL, leg.at=NULL, minmax=NULL, xlab="", ylab="", asp=x$asp, ...) {

	graphics::par(mar=x$mar)	
	
	plot(x$ext[1:2], x$ext[3:4], type = "n", xlab=xlab, ylab=ylab, asp=asp, ...)
	graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	
	if (x$leg$legend) {	
		if (inherits(leg.ext, "SpatExtent")) {
			leg.ext <- as.data.frame(rbind(as.vector(leg.ext)))
		} else if (is.numeric(leg.ext)) {
			leg.ext <- data.frame(rbind(leg.ext))
			if (ncol(leg.ext) != 4) {
				leg.ext <- NULL
			} else {
				colnames(leg.ext) <- c("xmin", "xmax", "ymin", "ymax")
			}
		}
		x$leg$ext <- leg.ext
		if (is.null(leg.levels)) {
			x$leg$levels <- 5
		} else {
			x$leg$levels <- leg.levels
		}
		x$leg$at <- leg.at
		x <- .plot.legend(x)
	}
	x
}	



.as.raster.continuous <- function(x, cols, minmax=NULL, ...) {
		
	out <- list()

	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA

	z <- stats::na.omit(as.vector(Z))
	if (length(z) == 0) stop("no values")
	if (is.null(minmax)) {
		minmax <- range(z)
	}
	interval <- (minmax[2]-minmax[1])/(length(cols)-1)
	breaks <- minmax[1] + interval * (0:(length(cols)-1))
		
	Z[] <- cols[as.integer(cut(Z, breaks, include.lowest=TRUE, right=FALSE))]
	out$raster <- as.raster(Z)

	out$leg <- list()
	out$leg$minmax <- minmax
	out$leg$cols <- cols
	out$leg$type <- "continuous"
	
	out
}


.prep.plot.data <- function(x, type, cols, maxcell, mar, leg, interpolate=FALSE, leg.shrink=c(0,0), leg.main="", leg.main.cex = 1, leg.digits, plot=FALSE, ...) {
	out <- list()
	out$mar <- mar
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
	} else {
		out$asp <- 1
	}
	out$ext <- as.vector(ext(x))

	if (type=="classes") {
		#ras <- .as.raster.classes(x, cols, ...)
	} else if (type=="range") {
		#ras <- .as.raster.range(x, cols, ...)
	} else {
		ras <- .as.raster.continuous(x, cols, ...)
	}
	out$r <- ras$r
	out$leg <- ras$leg
	if (is.na(leg) || isFALSE(leg)) {
		out$leg$legend <- FALSE
	} else {
		out$leg$legend <- TRUE
		out$leg$loc <- leg	

		if (missing(leg.digits)) {
			dif <- diff(out$leg$minmax)
			if (dif == 0) {
				leg.digits = 0;
			} else {
				leg.digits <- max(0, -floor(log10(dif/10)))
			}
		}
		out$leg$digits <- leg.digits
		out$leg$shrink <- leg.shrink
	}
	out$leg$main <- leg.main
	out$leg$main.cex <-  leg.main.cex	
	out$interpolate <- interpolate

	if (plot) {
		out <- .plotit(out, ...)
	}
	out
}



#setMethod("plot", signature(x="SpatRaster", y="numeric"), 

.plt <- function(x, y=1, col, type="continuous", mar=c(5.1, 4.1, 4.1, 7.1), legend="right", interpolate=FALSE, maxcell=50000, ...) {

		if (missing(col)) col <- rev(grDevices::terrain.colors(255))
		x <- x[[y]]
		if (!hasValues(x)) {
			stop("SpatRaster has no cell values")
		}
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		x <- .prep.plot.data(object, type=type, cols=col, mar=mar, leg=legend, interpolate, plot=TRUE, ...)
		invisible(x)
		#y <- .prep.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), leg="bottom")
		#y$interpolate <- F
	}
#}

#r <- rast(system.file("ex/test.tif", package="terra"))
#e <- c(177963, 179702, 333502, 333650) 
#.plt(r, leg="top", mar=c(2,2,2,2), leg.ext=e, leg.levels=3, leg.at=c(10, 666,1222), minmax=c(0,2000))
