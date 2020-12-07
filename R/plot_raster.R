

.as.raster.continuous <- function(out, x, type) {
		
	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	
	z <- stats::na.omit(as.vector(Z))
	n <- length(z)
	if (n == 0) stop("no values")
	if (type == "depends") {
		if (length(unique(z)) < 6) {
			return (.as.raster.classes(out, x))
		}
	} else if (length(unique(z)) == 1) {
		return (.as.raster.classes(out, x))
	}

	if (is.null(out$range)) {
		out$range <- range(z)
	} else {
		stopifnot(length(out$range) == 2)
		stopifnot(out$range[2] > out$range[1])
	}

	# avoid missing extreme values due to precision problems
	# this did not do it
	# Z <- round(Z, 10) 
	# so let's try this
	out$range[1] <- 0.9999999999 * out$range[1]
	out$range[2] <- 1.0000000001 * out$range[2]
	interval <- (out$range[2]-out$range[1])/(length(out$cols)-1)
	breaks <- out$range[1] + interval * (0:(length(out$cols)-1))
		
	Z[] <- out$cols[as.integer(cut(Z, breaks, include.lowest=TRUE, right=FALSE))]
	out$r <- as.raster(Z)

	out$legend_type <- "continuous"

	if (is.null(out$levels)) {
		out$levels <- 5
	} 
	if (is.null(out$leg$digits)) {
		dif <- diff(out$range)
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

	fz <- as.factor(Z)
	levs <- levels(fz)
	if (!is.null(out$levels)) {
		if (is.null(out$leg$legend)) {
			out$leg$legend <- as.character(out$levels)
		}
		levs <- out$levels
	} else if (!is.null(out$facts)) {
		uz <- as.integer(levels(fz))
		labs <- out$facts[[1]]$labels[uz+1]
		labs[is.na(labs)] <- "?"
		out$leg$legend <- labs
		out$levels <- uz
	} else {
		out$levels <- as.numeric(levs)
		if (is.null(out$leg$legend)) {
			out$leg$legend <- levs
		}
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
	
	if (is.null(out$leg$x)) {
		if (is.null(out$leg$ext)) {
			out$leg$x = "top"
			out$leg$y = NULL
		} else {
			if (length(out$leg$ext) == 4) {
				out$leg$x = out$leg$ext[1]
				out$leg$y = out$leg$ext[4]
			} else {
				out$leg$x = "top"
				out$leg$y = NULL
			}
		}	
	}	
	out
}


.as.raster.interval <- function(out, x, ...) {

	if (is.null(out$breaks)) {
		out$breaks <- 5
	} 

	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	fz <- cut(Z, out$breaks, include.lowest=TRUE, right=FALSE)

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


.as.raster.colortable <- function(out, x, ...) {
	z <- round(values(x))
	z[z<0 | z>255] <- NA
	z[is.nan(z) | is.infinite(z)] <- NA
	if (all(is.na(z))) {
		stop("no values")
	}
	out$cols <- grDevices::rgb(out$coltab[,1], out$coltab[,2], out$coltab[,3], out$coltab[,4], maxColorValue=255)
	z <- out$cols[z]
	z <- matrix(z, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
	out$r <- as.raster(z)
	out$legend_draw	 <- FALSE
	out
}


# legend 	
#	border="black", box.lwd = graphics::par("lwd"), box.lty = graphics::par("lty"), 
#	box.col = graphics::par("fg"), bty = "o", bg = graphics::par("bg"), xjust = 0, 
#	yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), text.width = NULL, 
#	text.col = graphics::par("col"), text.font = NULL, ncol = 1, horiz = FALSE, title = NULL,
 #   inset = 0, title.col = text.col, title.adj = 0.5, 

.plotit <- function(x, xlab="", ylab="", type = "n", yaxs="i", xaxs="i", asp=x$asp, axes=TRUE, new=NA, ...) {
	
	if ((!x$add) & (!x$legend_only)) {
	
		if (!any(is.na(x$mar))) { graphics::par(mar=x$mar) }
		plot(x$lim[1:2], x$lim[3:4], type=type, xlab=xlab, ylab=ylab, asp=asp, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
	}
	
	if (!x$legend_only) {
		graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
			angle = 0, interpolate = x$interpolate)	

		if (axes) x <- .plot.axes(x)			
	}
	if (x$legend_draw) {	
		if (x$legend_type == "continuous") {
			x <- do.call(.plot.cont.legend, list(x=x))
#		} else if (x$legend_type == "classes") {
		} else {
			#y <- do.call(.plot.class.legend, x$leg)
			x <- do.call(.plot.class.legend, x$leg)
		}
	}
	x
}	



.prep.plot.data <- function(x, type, maxcell, cols, mar, draw=FALSE, interpolate=FALSE,  
legend=TRUE, legend.only=FALSE, pax=list(), plg=list(), levels=NULL, add=FALSE,
 range=NULL, new=NA, breaks=NULL, coltab=NULL, facts=NULL, xlim=NULL, ylim=NULL, colNA=NA, alpha=NULL, ...) {

#mar=c(5.1, 4.1, 4.1, 7.1); legend=TRUE; axes=TRUE; pal=list(); pax=list(); maxcell=50000; draw=FALSE; interpolate=FALSE; legend=TRUE; legend.only=FALSE; pax=list(); pal=list(); levels=NULL; add=FALSE; range=NULL; new=NA; breaks=NULL; coltab=NULL; facts=NULL; xlim=NULL; ylim=NULL;
 
	out <- list()

	if (!(is.null(xlim) & is.null(ylim))) {
		e <- as.vector(ext(x))
		if (!is.null(xlim)) e[1:2] <- xlim
		if (!is.null(ylim)) e[3:4] <- ylim
		x <- crop(x, ext(e))
		out$ext <- as.vector(ext(x))
		out$lim <- e
	} else {
		out$lim <- out$ext <- as.vector(ext(x))
	}

	x <- spatSample(x, maxcell, method="regular", as.raster=TRUE)

	out$add <- isTRUE(add)
	out$mar <- mar

	out$axs <- pax 
	out$leg <- plg
	out$asp <- 1
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(out$ext[3:4]) * pi)/180)
	}
	
	if (!is.null(alpha)) {
		alpha <- clamp(alpha[1]*255, 0, 255)
		cols <- grDevices::rgb(t(grDevices::col2rgb(cols)), alpha=alpha, maxColorValue=255)
	} else {
		alpha <- 255
	}
	
	out$cols <- cols
	out$facts <- facts
	out$breaks <- breaks
	out$interpolate <- FALSE
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)

	if (type=="colortable") {
		out$coltab <- coltab
		out <- .as.raster.colortable(out, x)
	} else if (type=="classes") {
		out$levels <- levels
		out <- .as.raster.classes(out, x)
	} else if (type=="interval") {
		out <- .as.raster.interval(out, x)
	} else {
		out$interpolate <- isTRUE(interpolate)
		out$range <- range
		out <- .as.raster.continuous(out, x, type)
	}

	if (!is.null(colNA)) {
		if (!is.na(colNA)) {
			out$colNA <- grDevices::rgb(t(grDevices::col2rgb(colNA)), alpha=alpha, maxColorValue=255)
			out$r[is.na(out$r)] <- out$colNA
		}
	}

	if (draw) {
		out <- .plotit(out, new=new, ...)
	}
	invisible(out)
}


setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function(x, y=1, col, type, mar=c(5.1, 4.1, 4.1, 7.1), legend=TRUE, axes=TRUE, plg=list(), pax=list(), maxcell=50000, smooth=FALSE, range=NULL, levels=NULL, fun=NULL, colNA=NULL, alpha=NULL, ...) {

		x <- x[[y]]
		if (!hasValues(x)) { stop("SpatRaster has no cell values") }

		breaks <- list(...)$breaks
		coltab <- NULL
		facts  <- NULL
		if (!is.null(breaks)) {
			type <- "interval"
		} else {
			if (missing(type)) {
				if (x@ptr$hasColors()) {
					coltab <- cols(x)[[1]]
					type <- "colortable"
				} else if (is.factor(x)) {
					type <- "classes"
					facts <- levels(x)
				} else {
					type <- "depends"
				}
			} else {
				type <- match.arg(type, c("continuous", "classes", "interval"))
			}
		}
		
		if (missing(col)) col <- rev(grDevices::terrain.colors(25))
		x <- .prep.plot.data(x, type=type, maxcell=maxcell, cols=col, mar=mar, draw=TRUE, plg=plg, pax=pax, legend=isTRUE(legend), axes=isTRUE(axes), coltab=coltab, facts=facts, interpolate=smooth, levels=levels, range=range, ...)

		if (!is.null(fun)) {
			if (!is.null(formals(fun))) {
				fun(1)
			} else {
				fun()
			}
		}
		invisible(x)
	}
)

setMethod("plot", signature(x="SpatRaster", y="missing"), 
	function(x, y, maxcell=50000, main, nc, nr, maxnl=16, ...)  {

		nl <- max(1, min(nlyr(x), maxnl))

		if (nl==1) {
			if (missing(main)) {
				out <- plot(x, 1, maxcell=maxcell, ...)
			} else {
				out <- plot(x, 1, maxcell=maxcell, main=main[1], ...)
			}
			return(invisible(out))
		}

		nrnc <- .get_nrnc(nr, nc, nl)
		old.par <- graphics::par(no.readonly = TRUE) 
		on.exit(graphics::par(old.par))
		graphics::par(mfrow=nrnc, mar=c(2, 2, 2, 4))
		maxcell=maxcell/(nl/2)
			
		if (missing("main")) {
			main <- names(x)
		} else {
			main <- rep_len(main, nl)	
		}
		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		for (i in 1:nl) {
			plot(x, i, main=main[i], ...)
		}
	}
)





#mar=c(5.1, 4.1, 4.1, 7.1); legend=TRUE; axes=TRUE; plg=list(); pax=list(); maxcell=50000

#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
#x <- .prep.plot.data(object, type="classes", cols=rainbow(25), mar=rep(3,4), draw=T)
#x <- .prep.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), draw=T)

#r <- rast(system.file("ex/test.tif", package="terra"))
#plot(r)
#e <- c(177963, 179702, 333502, 333650) 
#plot(r, mar=c(3,3,3,3), plg=list(loc="top", ext=e, levels=3, at=c(10, 666,1222), range=c(0,2000)))
#plot(r, type="interval", levels=c(0,500,3000), plg=list(legend=c("low", "high"), x="topleft"))
#plot(r, type="interval", plg=list(x="topleft"))
#plot(r, type="interval", plg=list(cex=.8, bty="n"), pax=list(cex.axis=.8, las=1))
#plot(r, type="interval", plg=list(cex=.8, bty="n", ncol=2, x=178000, y=335250), pax=list(cex.axis=.8, las=1))
 
#par(mfrow=c(1,2))
#plot(r, type="interval", mar=c(2,4,2,0), plg=list(inset=0.05, x="topleft"), pax=list(sides=c(1,2), cex.axis=0.8))
#plot(r, type="interval", mar=c(2,0,2,4), legend=FALSE, pax=list(sides=c(3,4), cex.axis=0.8))

#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
#type="classes"; cols=rainbow(25); mar=rep(3,4); draw=TRUE
 
#r <- rast(system.file("ex/test.tif", package="terra"))
#plot(r, type="interval", mar=c(2,4,2,0), plg=list(x="topleft", inset=0.05), pax=list(sides=c(1,2), cex.axis=0.8))
#plot(r)

# f <- system.file("ex/test.tif", package="terra") 
# r <- rast(f)
# plot(r)
# plot(r, type="interval")
# plot(r, type="classes")

# d <- (r > 400) + (r > 600)
# plot(d)
# e <- c(178000,178200,332000,333000)
# plot(d, plg=list(ext=e, title="Title\n", title.cex=2))

# plot(d, type="interval", levels=0:3) 
# plot(d, type="interval", levels=3, plg=list(legend=c("0-1", "1-2", "2-3"))) 
# plot(d, type="classes", plg=list(legend=c("M", "X", "A")))

# r <- rast(f)
# x <- trunc(r/600)
# x <- as.factor(x)
# levels(x) <- c("earth", "wind", "fire")
# plot(x)

