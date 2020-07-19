

.plotit <- function(x, minmax=NULL, xlab="", ylab="", type = "n", asp=x$asp, ...) {
	
	if (!is.na(x$mar)) graphics::par(mar=x$mar)	
	
	plot(x$ext[1:2], x$ext[3:4], type=type, xlab=xlab, ylab=ylab, asp=asp, ...)
	
	graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	
	
	if (x$leg$draw) {	
		x <- .plot.legend(x)
	}
	x
}	



.as.raster.continuous <- function(out, x, minmax=NULL, ...) {
		
	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA

	z <- stats::na.omit(as.vector(Z))
	if (length(z) == 0) stop("no values")
	if (is.null(minmax)) {
		minmax <- range(z)
	}
	interval <- (minmax[2]-minmax[1])/(length(out$cols)-1)
	breaks <- minmax[1] + interval * (0:(length(out$cols)-1))
		
	Z[] <- out$cols[as.integer(cut(Z, breaks, include.lowest=TRUE, right=FALSE))]
	out$r <- as.raster(Z)

	out$leg$minmax <- minmax
	out$leg$type <- "continuous"

	if (is.null(out$leg$levels)) {
		out$leg$levels <- 5
	} 
	if (is.null(out$leg$digits)) {
		dif <- diff(out$leg$minmax)
		if (dif == 0) {
			out$leg$digits = 0;
		} else {
			out$leg$digits <- max(0, -floor(log10(dif/10)))
		}
	}
	
	if (is.null(out$leg$loc)) out$leg$loc <- "right"
	out
}


.as.raster.classes <- function(out, x, ...) {

	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	fz <- as.factor(Z)
	nlevs <- length(levels(fz))
	if (nlevs == 0) {
		stop("no values")
	}
	cols <- rep_len(cols, nlevs)
	Z[] <- cols[as.numeric(fz)]
	
	out$r <- as.raster(Z)

	out$leg$cols <- cols
	out$leg$levels <- as.numeric(levels(fz))
	out$leg$labels <- levels(fz)
	
	out$leg$type <- "classes"
	out
}


.prep.plot.data <- function(x, type, cols, mar, draw=FALSE, interpolate=FALSE,  
legend=TRUE, leg.shrink=c(0,0), leg.main="", leg.main.cex = 1, leg.digits=NULL, leg.loc=NULL, leg.ext=NULL, leg.levels=NULL, leg.at=NULL, ...) {

	out <- list()
	out$mar <- mar
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
	} else {
		out$asp <- 1
	}
	out$ext <- as.vector(ext(x))
	out$cols <- cols
	out$interpolate <- isTRUE(interpolate)

	out$leg$loc <- leg.loc
	out$leg$digits <- leg.digits
	out$leg$draw <- isTRUE(legend)
	out$leg$shrink <- leg.shrink
	out$leg$main <- leg.main
	out$leg$main.cex <- leg.main.cex	
	out$leg$at <- leg.at
	out$leg$ext <- as.vector(leg.ext)

	if (type=="classes") {
		out <- .as.raster.classes(out, x, ...)
	} else if (type=="range") {
		#ras <- .as.raster.range(x, cols, ...)
	} else {
		out <- .as.raster.continuous(out, x, ...)
	}

	if (draw) {
		out <- .plotit(out, ...)
	}
	out
}



#setMethod("plot", signature(x="SpatRaster", y="numeric"), 

.plt <- function(x, y=1, col, type="continuous", mar=c(5.1, 4.1, 4.1, 7.1), maxcell=50000, ...) {

		if (missing(col)) col <- rev(grDevices::terrain.colors(255))
		x <- x[[y]]
		if (!hasValues(x)) {
			stop("SpatRaster has no cell values")
		}
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		x <- .prep.plot.data(object, type=type, cols=col, mar=mar, draw=TRUE, ...)
		invisible(x)
		#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
		#x <- .prep.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), draw=T)
	}
#}

#r <- rast(system.file("ex/test.tif", package="terra"))
#e <- c(177963, 179702, 333502, 333650) 
#terra:::.plt(r, leg.loc="top", mar=c(2,2,2,2), leg.ext=e, leg.levels=3, leg.at=c(10, 666,1222), minmax=c(0,2000))
