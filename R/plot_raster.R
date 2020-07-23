
.plot.axes <- function(x) {
	if (!is.null(x$axs$sides)) {
		usr <- graphics::par("usr")
		sides <- x$axs$sides
		x$axs$sides <- NULL
		sides <- round(unique(sides))
		sides[sides > 1 & sides < 5]
		for (s in sides) {
			if (s %in% c(1,3)) {
				ur <- usr[2] - usr[1]
				at <- c(usr[1]-10*ur, usr[2]+10*ur)
			} else {
				ur <- usr[4] - usr[3]
				at <- c(usr[3]-10*ur, usr[4]+10*ur)
			}
			graphics::axis(s, at=at, labels=c("",""), lwd.ticks=0)
			x$axs$side <- s
			do.call(graphics::axis, x$axs)
		}
		x$axs$sides <- x$sides
	} else {
		x$axs$side <- 1
		do.call(graphics::axis, x$axs)
		x$axs$side <- 2
		do.call(graphics::axis, x$axs)
		graphics::box()
	}
	x$axs$side <- NULL
	x
}



# legend 	
#	border="black", box.lwd = graphics::par("lwd"), box.lty = graphics::par("lty"), 
#	box.col = graphics::par("fg"), bty = "o", bg = graphics::par("bg"), xjust = 0, 
#	yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), text.width = NULL, 
#	text.col = graphics::par("col"), text.font = NULL, ncol = 1, horiz = FALSE, title = NULL,
 #   inset = 0, title.col = text.col, title.adj = 0.5, 

.plotit <- function(x, xlab="", ylab="", type = "n", asp=x$asp, axes=TRUE, ...) {
	
	if (!x$add) {
		if (!any(is.na(x$mar))) { graphics::par(mar=x$mar) }
		plot(x$ext[1:2], x$ext[3:4], type=type, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, ...)
		if (axes) .plot.axes(x)
	}
	
	graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	
	
	if (x$legend_draw) {	
		if (x$legend_type == "continuous") {
			x <- do.call(.plot.cont.legend, list(x=x))
#		} else if (x$legend_type == "classes") {
		} else {
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
		if (is.null(out$leg$legend)) {
			out$leg$legend <- levs
		}
	} else {
		if (is.null(out$leg$legend)) {
			out$leg$legend <- as.character(out$levels)
		}
		levs <- out$levels
	}
	stopifnot(length(out$leg$legend) == length(out$levels))
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
	out$legend_type <- "classes"
	
	#&& is.null(out$leg$ext)) 
	if (is.null(out$leg$x)) {
		out$leg$x = "top"
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
		out$leg$x <- "top"
	}
	out
}

# leg.shrink=c(0,0), leg.main=NULL, leg.main.cex = 1, leg.digits=NULL, leg.loc=NULL, leg.ext=NULL, leg.levels=NULL, leg.labels=NULL, leg.at=NULL, 

.prep.plot.data <- function(x, type, cols, mar, draw=FALSE, interpolate=FALSE,  
legend=TRUE, pax=list(), pal=list(), levels=NULL, add=FALSE, ...) {

	out <- list()
	out$add <- isTRUE(add)
	out$mar <- mar
	out$ext <- as.vector(ext(x))
	out$axs <- pax 
	out$leg <- pal
	out$asp <- 1
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(out$ext[3:4]) * pi)/180)
	}
	out$cols <- cols
	out$levels <- levels
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
		out <- .plotit(out, ...)
	}
	out
}


setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function(x, y=1, col, type=c("continuous", "classes", "interval"), mar=c(5.1, 4.1, 4.1, 7.1), legend=TRUE, axes=TRUE, pal=list(), pax=list(), maxcell=50000, ...) {

		type <- match.arg(type)
		if (!hasValues(x)) { stop("SpatRaster has no cell values") }
		x <- x[[y]]
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		if (missing(col)) col <- rev(grDevices::terrain.colors(25))
		x <- .prep.plot.data(object, type=type, cols=col, mar=mar, draw=TRUE, pal=pal, pax=pax, legend=isTRUE(legend), axes=isTRUE(axes), ...)
		invisible(x)
	}
)


#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
#x <- .prep.plot.data(object, type="classes", cols=rainbow(25), mar=rep(3,4), draw=T)
#x <- .prep.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), draw=T)

#r <- rast(system.file("ex/test.tif", package="terra"))
#map(r)
#e <- c(177963, 179702, 333502, 333650) 
#map(r, mar=c(3,3,3,3), pal=list(loc="top", ext=e, levels=3, at=c(10, 666,1222), range=c(0,2000)))
#map(r, type="interval", levels=c(0,500,3000), pal=list(legend=c("low", "high"), x="topleft"))
#map(r, type="interval", pal=list(x="topleft"))
#map(r, type="interval", pal=list(cex=.8, bty="n"), pax=list(cex.axis=.8, las=1))
#map(r, type="interval", pal=list(cex=.8, bty="n", ncol=2, x=178000, y=335250), pax=list(cex.axis=.8, las=1))
 
#par(mfrow=c(1,2))
#map(r, type="interval", mar=c(2,4,2,0), pal=list(inset=0.05, x="topleft"), pax=list(sides=c(1,2), cex.axis=0.8))
#map(r, type="interval", mar=c(2,0,2,4), legend=FALSE, pax=list(sides=c(3,4), cex.axis=0.8))

#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
#type="classes"; cols=rainbow(25); mar=rep(3,4); draw=TRUE
 
#r <- rast(system.file("ex/test.tif", package="terra"))
#map(r, type="interval", mar=c(2,4,2,0), pal=list(x="topleft", inset=0.05), pax=list(sides=c(1,2), cex.axis=0.8))
#map(r)

# f <- system.file("ex/test.tif", package="terra") 
# r <- rast(f)
# map(r)
# map(r, type="interval")
# map(r, type="classes")

# d <- (r > 400) + (r > 600)
# map(d)
# e <- c(178000,178200,332000,333000)
# map(d, pal=list(ext=e, title="Title\n", title.cex=2))

# map(d, type="interval", levels=0:3) 
# map(d, type="interval", levels=3, pal=list(legend=c("0-1", "1-2", "2-3"))) 
# map(d, type="classes", pal=list(legend=c("M", "X", "A")))

# r <- rast(f)
# x <- trunc(r/600)
# x <- as.factor(x)
# levels(x) <- c("earth", "wind", "fire")
# plot(x)
