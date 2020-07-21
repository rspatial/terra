
.plot.axes <- function(x) {
	if (!is.null(x$axs$axes)) {
		usr <- graphics::par("usr")
		axarg <- x$axs$axes
		x$axs$axes <- NULL
		axarg <- round(unique(axarg))
		axarg[axarg > 1 & axarg < 5]
		for (i in axarg) {
			if (i %in% c(1,3)) {
				ur <- usr[2] - usr[1]
				at = c(usr[1]-10*ur, usr[2]+10*ur)
			} else {
				ur <- usr[4] - usr[3]
				at = c(usr[3]-10*ur, usr[4]+10*ur)
			}
			graphics::axis(i, at=at, labels=c("",""), lwd.ticks=0)
			x$axs$side <- i
			do.call(graphics::axis, x$axs)
			#graphics::axis(i, cex.axis=cex.axis, las=las, ...)
		}
		x$axs <- x$axarg
	} else {
		graphics::axis(1)
		graphics::axis(2)
		graphics::box()
	}
	x
}



# legend 	
#	border="black", box.lwd = graphics::par("lwd"), box.lty = graphics::par("lty"), 
#	box.col = graphics::par("fg"), bty = "o", bg = graphics::par("bg"), xjust = 0, 
#	yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), text.width = NULL, 
#	text.col = graphics::par("col"), text.font = NULL, ncol = 1, horiz = FALSE, title = NULL,
 #   inset = 0, title.col = text.col, title.adj = 0.5, 

.plotit <- function(x, xlab="", ylab="", type = "n", asp=x$asp, axes, ...) {
	
	if (!any(is.na(x$mar))) {
		graphics::par(mar=x$mar)	
	}
	
	plot(x$ext[1:2], x$ext[3:4], type=type, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, ...)

	graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	

	.plot.axes(x)
	
	if (x$legend_draw) {	
	
		if (x$legend_type == "continuous") {
			if (is.null(x$legend_ext)) {
				x <- .get.leg.extent(x)
			} else {
				x <- .get.leg.coords(x)	
			}
			x <- do.call(.plot.cont.legend, list(x=x))
			
		} else if (x$legend_type == "classes") {
			y <- do.call(.plot.class.legend, x$leg)
		}
	}
	x
}	



.as.raster.continuous <- function(out, x, ...) {
		
	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA

	z <- stats::na.omit(as.vector(Z))
	if (length(z) == 0) stop("no values")
	if (is.null(out$leg$range)) {
		out$leg$range <- range(z)
	}
	interval <- (out$leg$range[2]-out$leg$range[1])/(length(out$cols)-1)
	breaks <- out$leg$range[1] + interval * (0:(length(out$cols)-1))
		
	Z[] <- out$cols[as.integer(cut(Z, breaks, include.lowest=TRUE, right=FALSE))]
	out$r <- as.raster(Z)

	out$legend_type <- "continuous"

	if (is.null(out$levels)) {
		out$levels <- 5
	} 
	if (is.null(out$leg$digits)) {
		dif <- diff(out$leg$range)
		if (dif == 0) {
			out$leg_digits = 0;
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
	if (all(is.na(Z))) {
		stop("no values")
	}

	if (is.null(out$levels)) {
		fz <- as.factor(Z)
		levs <- levels(fz)
		out$levels <- as.numeric(levs)
		out$leg$legend <- levs
	} else {
		if (is.null(out$leg$legend)) {
			out$leg$legend <- as.character(out$levels)
		} else {
			stopifnot(length(out$leg$legend) == length(out$levels))
		}
		levs <- out$levels
	}
	nlevs <- length(levs)
	
	cols <- out$cols
	ncols <- length(cols)
	if (nlevs < ncols) {
		i <- trunc((ncols / nlevs) * 1:nlevs)
		cols <- cols[i]
	} else {
		cols <- rep_len(cols, nlevs)
	}
	out$leg$fill <- cols
	Z[] <- cols[as.numeric(fz)]
	out$r <- as.raster(Z)
	out$leg$type <- "classes"
	
	#&& is.null(out$leg$ext)) 
	if (is.null(out$leg$x)) {
		out$leg$x <- "topright"
	}
	out
}


.as.raster.interval <- function(out, x, ...) {

	if (is.null(out$levels)) {
		out$levels <- 5
	} 

	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	fz <- cut(Z, out$levels, include.lowest=TRUE, right=FALSE)

	levs <- levels(fz)
	nlevs <- length(levs)
	
	cols <- out$cols
	ncols <- length(cols)
	if (nlevs < ncols) {
		i <- trunc((ncols / nlevs) * 1:nlevs)
		cols <- cols[i]
	} else {
		cols <- rep_len(cols, nlevs)
	}
	#out$cols <- cols
	out$leg$fill <- cols
	Z[] <- cols[as.integer(fz)]
	out$r <- as.raster(Z)
	#out$leg$levels <- levels(fz)
	out$legend_type <- "classes"
	
	if (!is.null(out$leg$legend)) {
		stopifnot(length(out$leg$legend) == nlevs)	
	} else {
		levs <- gsub("]", "", gsub(")", "", gsub("\\[", "", levs)))
		levs <- paste(levs, collapse=",")
		m <- matrix(as.numeric(unlist(strsplit(levs, ","))), ncol=2, byrow=TRUE)
		m <- apply(m, 1, function(i) paste(i, collapse=" - "))
		out$leg$legend <- m
	}
	
	if (is.null(out$leg$x)) { # && is.null(out$leg$ext)) {
		out$leg$x <- "topright"
	}
	out
}

# leg.shrink=c(0,0), leg.main=NULL, leg.main.cex = 1, leg.digits=NULL, leg.loc=NULL, leg.ext=NULL, leg.levels=NULL, leg.labels=NULL, leg.at=NULL, 

.prep.plot.data <- function(x, type, cols, mar, draw=FALSE, interpolate=FALSE,  
legend=TRUE, pax=list(), plg=list(), ...) {

	out <- list()
	out$mar <- mar
	out$axs <- pax 
	out$leg <- plg
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
	} else {
		out$asp <- 1
	}
	out$ext <- as.vector(ext(x))
	out$cols <- cols
	out$interpolate <- isTRUE(interpolate)
	out$legend_draw <- isTRUE(legend)

	if (type=="classes") {
		out <- .as.raster.classes(out, x)
	} else if (type=="interval") {
		out <- .as.raster.interval(out, x)
	} else {
		out <- .as.raster.continuous(out, x)
	}

	if (draw) {
		out <- .plotit(out, pax, ...)
	}
	out
}



#setMethod("plot", signature(x="SpatRaster", y="numeric"), 

map <- function(x, y=1, col, type="continuous", mar=c(5.1, 4.1, 4.1, 7.1), maxcell=50000, plg=list(), pax=list(), ...) {

		if (missing(col)) col <- rev(grDevices::terrain.colors(25))
		x <- x[[y]]
		if (!hasValues(x)) {
			stop("SpatRaster has no cell values")
		}
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		x <- .prep.plot.data(object, type=type, cols=col, mar=mar, draw=TRUE, plg=plg, pax=pax, ...)
		invisible(x)
		#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
		#x <- .prep.plot.data(object, type="classes", cols=rainbow(25), mar=rep(3,4), draw=T)
		#x <- .prep.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), draw=T)
	}
#}

#r <- rast(system.file("ex/test.tif", package="terra"))
#e <- c(177963, 179702, 333502, 333650) 
#map(r)

#map(r, mar=c(2,2,4,3), plg=list(loc="top", ext=e, levels=3, at=c(10, 666,1222), range=c(0,2000)))

#.plt(r, type="interval", leg.lev=c(0,500,3000), leg.lab=c("low", "high"), leg.loc="topleft")

#map(r, type="interval", plg=list(x="topleft"))

#par(mfrow=c(1,2))
#map(r, type="interval", mar=c(2,4,2,0), plg=list(inset=0.05, x="topleft"), pax=list(axes=c(1,2), cex.axis=0.8))
#map(r, type="interval", mar=c(2,0,2,4), legend=FALSE, pax=list(axes=c(3,4), cex.axis=0.8))
#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
#type="classes"; cols=rainbow(25); mar=rep(3,4); draw=TRUE
 
#r <- rast(system.file("ex/test.tif", package="terra"))
#map(r, type="interval", mar=c(2,4,2,0), plg=list(x="topleft", inset=0.05), pax=list(axes=c(1,2), cex.axis=0.8))
#map(r)


