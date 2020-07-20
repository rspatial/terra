

.plotit <- function(x, minmax=NULL, xlab="", ylab="", type = "n", asp=x$asp, 
	cex=1, cex.axis=1,  cex.lab=1, main="", cex.main=1, sub="", cex.sub=1, las=0,

# legend 	
	border="black", box.lwd = graphics::par("lwd"), box.lty = graphics::par("lty"), 
	box.col = graphics::par("fg"), bty = "o", bg = graphics::par("bg"), xjust = 0, 
	yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), text.width = NULL, 
	text.col = graphics::par("col"), text.font = NULL, ncol = 1, horiz = FALSE, title = NULL,
    inset = 0, title.col = text.col, title.adj = 0.5, 
	
	...) {
	
	if (!any(is.na(x$mar))) {
		graphics::par(mar=x$mar)	
	}
	
	axarg <- NULL
	if (is.logical(x$axes)) {
		axes <- x$axes
	} else {
		axarg <- x$axes
		axes <- FALSE
	}
	plot(x$ext[1:2], x$ext[3:4], type=type, xlab=xlab, ylab=ylab, asp=asp, cex.axis=cex.axis, axes=axes, main=main, cex=cex, cex.main=cex.main, sub=sub, cex.sub=cex.sub, ...)

	if (!is.null(axarg)) {
		usr <- graphics::par("usr")
		axarg <- round(unique(axarg))
		axarg[axarg > 1 & axarg < 5]
		for (i in axarg) {
			graphics::axis(i, cex.axis=cex.axis, las=las)
			if (i %in% c(1,3)) {
				ur <- usr[2] - usr[1]
				at = c(usr[1]-10*ur, usr[2]+10*ur)
			} else {
				ur <- usr[4] - usr[3]
				at = c(usr[3]-10*ur, usr[4]+10*ur)
			}
			graphics::axis(i, at=at, labels=c("",""), lwd.ticks=0)
		}
	}
	
	graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], 
		angle = 0, interpolate = x$interpolate)	
	
	if (x$leg$draw) {	
	
		if (x$leg$type == "continuous") {
			if (is.null(x$leg$ext)) {
				x <- .get.leg.extent(x)
			} else {
				x <- .get.leg.coords(x)	
			}
			x <- .plot.cont.legend(x, ...)
		} else if (x$leg$type == "classes") {
			x <- .plot.class.legend(x, cex=cex, border=border, box.lwd=box.lwd, box.lty=box.lty, box.col=box.col, bty=bty, bg=bg, xjust=xjust, yjust=yjust, x.intersp=x.intersp, y.intersp=y.intersp, adj=adj, text.width=text.width, text.col=text.col, text.font=text.font, ncol=ncol, horiz=horiz, title=title, inset=inset, title.col=title.col, title.adj=title.adj)
		}
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
	if (all(is.na(Z))) {
		stop("no values")
	}

	if (is.null(out$leg$levels)) {
		fz <- as.factor(Z)
		levs <- levels(fz)
		out$leg$levels <- as.numeric(levs)
		out$leg$labels <- levs
	} else {
		if (is.null(out$leg$lables)) {
			out$leg$labels <- as.character(out$leg$levels)
		} else {
			stopifnot(length(out$leg$lables) == length(out$leg$levels))
		}
		levs <- out$leg$levels
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
	out$cols <- cols
	Z[] <- cols[as.numeric(fz)]
	out$r <- as.raster(Z)
	out$leg$type <- "classes"
	
	if (is.null(out$leg$loc) && is.null(out$leg$ext)) {
		out$leg$loc <- "topright"
	}
	out
}


.as.raster.range <- function(out, x, ...) {

	if (is.null(out$leg$levels)) {
		out$leg$levels <- 5
	} 

	Z <- as.matrix(x, TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	fz <- cut(Z, out$leg$levels, include.lowest=TRUE, right=FALSE)

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
	out$cols <- cols
	Z[] <- cols[as.integer(fz)]
	out$r <- as.raster(Z)
	#out$leg$levels <- levels(fz)
	out$leg$type <- "classes"
	
	if (!is.null(out$leg$labels)) {
		stopifnot(length(out$leg$labels) == nlevs)	
	} else {
		levs <- gsub("]", "", gsub(")", "", gsub("\\[", "", levs)))
		levs <- paste(levs, collapse=",")
		m <- matrix(as.numeric(unlist(strsplit(levs, ","))), ncol=2, byrow=TRUE)
		m <- apply(m, 1, function(i) paste(i, collapse=" - "))
		out$leg$labels <- m
	}
	
	if (is.null(out$leg$loc) && is.null(out$leg$ext)) {
		out$leg$loc <- "topright"
	}
	out
}


.prep.plot.data <- function(x, type, cols, mar, axes=TRUE, draw=FALSE, interpolate=FALSE,  
legend=TRUE, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex = 1, leg.digits=NULL, leg.loc=NULL, leg.ext=NULL, leg.levels=NULL, leg.labels=NULL, leg.at=NULL, ...) {

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
	out$axes <- axes
	out$interpolate <- isTRUE(interpolate)

	out$leg$loc <- leg.loc
	out$leg$digits <- leg.digits
	out$leg$levels <- leg.levels
	out$leg$labels <- leg.labels
	out$leg$draw   <- isTRUE(legend)
	out$leg$shrink <- leg.shrink
	out$leg$main$text <- leg.main
	out$leg$main$cex  <- leg.main.cex	
	out$leg$at  <- leg.at
	out$leg$ext <- as.vector(leg.ext)

	if (type=="classes") {
		out <- .as.raster.classes(out, x, ...)
	} else if (type=="range") {
		out <- .as.raster.range(out, x, ...)
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

		if (missing(col)) col <- rev(grDevices::terrain.colors(25))
		x <- x[[y]]
		if (!hasValues(x)) {
			stop("SpatRaster has no cell values")
		}
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		x <- .prep.plot.data(object, type=type, cols=col, mar=mar, draw=TRUE, ...)
		invisible(x)
		#object <- spatSample(r, Inf, method="regular", as.raster=TRUE)
		#x <- .prep.plot.data(object, type="classes", cols=rainbow(25), mar=rep(3,4), draw=T)
		#x <- .prep.plot.data(object, type="continuous", cols=rainbow(25), mar=rep(3,4), draw=T)
	}
#}

#r <- rast(system.file("ex/test.tif", package="terra"))
#e <- c(177963, 179702, 333502, 333650) 
#.plt(r, type="range")
#.plt(r, leg.loc="top", mar=c(2,2,2,2), leg.ext=e, leg.levels=3, leg.at=c(10, 666,1222), minmax=c(0,2000))
#.plt(r, type="range", leg.lev=c(0,500,3000), leg.lab=c("low", "high"), leg.loc="topleft")

#par(mfrow=c(1,2))
#.plt(r, type="range", mar=c(2,4,2,0), leg.loc="topleft", axes=c(1,2), inset=0.05, cex.axis=0.8)
#box()
#.plt(r, type="range", mar=c(2,0,2,4), legend=FALSE, axes=c(3,4), cex.axis=0.8)
#box()
 
