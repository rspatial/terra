# Author: Robert J. Hijmans
# Date :  June 2019
# Version 1.0
# License GPL v3


.plotLines <- function(x, cols, ...) {
	g <- geom(x)
	g <- split(g, g[,1])
	g <- lapply(g, function(x) split(x, x[,2]))
	#p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4], ...))	
	for (i in 1:length(g)) {
		x <- g[[i]]
		for (j in 1:length(x)) {
			lines(x[[j]][,3:4], col=cols[i])
		}
	}
}

.plotPolygons <- function(x, cols, border=NULL, density=NULL, angle=45, ...) {
	g <- geom(x)
	g <- split(g, g[,1])
	g <- lapply(g, function(y) split(y, y[,2]))
	if (!is.null(border)) {
		border <- rep_len(border, length(g))
	} else {
		border <- NA
	}
	if (!is.null(density)) {
		density <- rep_len(density, length(g))
		angle <- rep_len(angle, length(g))
	}
	for (i in 1:length(g)) {
		gg <- g[[i]]
		for (j in 1:length(gg)) {
			a <- gg[[j]]
			if (any(a[,5] > 0)) {
				a <- split(a, a[,5]) 
				a <- lapply(a, function(i) rbind(i, NA))
				a <- do.call(rbind, a )
				a <- a[-nrow(a), ]
				# g[[i]][[1]] <- a 
			}
			if (!is.null(density)) {
				graphics::polypath(a[,3:4], col=NA, rule = "evenodd", border=border[i], ...)
				graphics::polygon(a[,3:4], col=cols[i], density=density[i], angle=angle[i], border=NA, ...)
			} else {
				graphics::polypath(a[,3:4], col=cols[i], rule = "evenodd", border=border[i], ...)
			}
		}
	}
}
.getCols <- function(n, cols) {
	if (!is.null(cols)) {
		ncols <- length(cols)
		if (ncols > n) {
			steps <- ncols/n
			i <- round(seq(1, ncols, steps))
			cols <- cols[i]				
		} else if (ncols < n) {
			cols <- rep_len(cols, n)
		}
	} 
	cols
}


.vect.legend.classes <- function(out) {
	ucols <- .getCols(length(out$uv), out$cols)
	uv <- sort(out$uv)
	i <- match(out$v, out$uv)
	out$cols <- ucols[i]

	out$levels <- out$uv
	out$leg$legend <- out$uv
	nlevs <- length(out$uv)
		
	cols <- out$cols
	ncols <- length(cols)
	if (nlevs < ncols) {
		i <- trunc((ncols / nlevs) * 1:nlevs)
		cols <- cols[i]
	} else {
		cols <- rep_len(cols, nlevs)
	}
	out$leg$fill <- cols
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


.vect.legend.continuous <- function(out) {


	z <- stats::na.omit(out$v)
	n <- length(z)
	if (n == 0) stop("no values")
	if (out$legend_type == "depends") {
		if (length(unique(z)) < 6) {
			return (.vect.legend.classes(out))
		}
	} else if (length(unique(z)) == 1) {
		return (.vect.legend.classes(out))
	}
	out$range <- range(z)
	
	interval <- (out$range[2]-out$range[1])/(length(out$cols)-1)
	breaks <- out$range[1] + interval * (0:(length(out$cols)-1))
		
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


.vect.legend.interval <- function(out, ...) {

	if (is.null(out$breaks)) {
		out$breaks <- 5
	} 
	fz <- cut(out$v, out$breaks, include.lowest=TRUE, right=FALSE)
	out$vcut <- as.integer(fz)
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
	out$leg$fill <- cols
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



.vplot <- function(x, col, add=FALSE, border="black", xlab="", ylab="", asp=NULL, density=NULL, angle=45, ...) {
	
	axes=FALSE 
	
	gtype <- geomtype(x)
	if (is.null(asp)) {
		if (isLonLat(x, perhaps=TRUE, warn=FALSE)) {
				asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} else {
			asp <- 1
		}
	}
	if (gtype == "points") {
		if (missing(col)) col = "black"
		g <- geom(x)
		if (add) {
			points(g[,3], g[,4], col=col, ...)			
		} else {
			plot(g[,3], g[,4], col=col, axes=axes, xlab=xlab, ylab=ylab, asp=asp, ...)
		}
	} else {
		e <- matrix(as.vector(ext(x)), 2)
		#if (!add) {
		#	plot(e, type="n", axes=axes, xlab=xlab, ylab=ylab, asp=asp, ...)
		#}
		if (gtype == "polygons") {
			.plotPolygons(x, col, border=border, density=density, angle=angle, ...)
		} else {
			if (is.null(col)) col = rep("black", size(x))
			.plotLines(x, col, ...)
		}
	}
}


.plot.vect.map <- function(x, out, xlab="", ylab="", type = "n", yaxs="i", xaxs="i", asp=out$asp, new=NA, ...) {
	
	if ((!out$add) & (!out$legend_only)) {
		if (!any(is.na(out$mar))) { graphics::par(mar=out$mar) }
		plot(out$lim[1:2], out$lim[3:4], type="n", xlab=xlab, ylab=ylab, asp=asp, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
	}

	if (out$legend_draw) {	
		if (out$legend_type == "classes") {
			out <- .vect.legend.classes(out)
			cols <- out$cols
		} else if (out$legend_type == "interval") {
			out <- .vect.legend.interval(out)
			cols <- out$cols[out$vcut]
		} else {
			out <- .vect.legend.continuous(out)
			brks <- seq(min(out$v, na.rm=TRUE), max(out$v, na.rm=TRUE), length.out = length(out$cols))
			grps <- cut(out$v, breaks = brks, include.lowest = TRUE)
			cols <- out$cols[grps]
		}				
		if (!out$legend_only) {
			.vplot(x, cols, add=out$add, ...) 
		}
	} else {
		.vplot(x, out$cols, add=out$add, ...) 
	}

	if (out$axes & !out$add) {
		out <- .plot.axes(out)	
	}

	if (out$legend_draw) {	
		if (out$legend_type == "continuous") {
			out <- do.call(.plot.cont.legend, list(x=out))
		} else {
			out <- do.call(.plot.class.legend, out$leg)
		}
	}
	out
}	



.prep.vect.data <- function(x, y, type="depends", cols, mar, draw=FALSE,   
legend=TRUE, legend.only=FALSE, pax=list(), plg=list(), levels=NULL, add=FALSE,
 range=NULL, new=NA, breaks=NULL, coltab=NULL, facts=NULL, xlim=NULL, ylim=NULL, colNA=NA, alpha=NULL, axes=TRUE, ...) {

	out <- list()

	if (!(is.null(xlim) & is.null(ylim))) {
		e <- as.vector(ext(x))
		if (!is.null(xlim)) e[1:2] <- xlim
		if (!is.null(ylim)) e[3:4] <- ylim
		out$ext <- as.vector(ext(x))
			out$lim <- e
	} else {
		out$lim <- out$ext <- as.vector(ext(x))
	}
	out$add <- isTRUE(add)
	out$axes <- isTRUE(axes)
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

# we may need to loop over y
	y = y[1]
	
	out$v <- unlist(x[, y, drop=TRUE], use.names=FALSE)
	out$uv <- unique(out$v)
	if (missing(type)) {
		if (!is.numeric(out$uv) | length(out$uv) < 10) {
			type <- "classes"
		} else {
			type <- "continuous"
		}
	} else {
		type <- match.arg(type, c("continuous", "classes", "interval", "depends", "none"))
	}
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)
	if (type=="none") {
		out$legend_draw <- FALSE
	} else if (type=="classes") {
		out$levels <- levels
	} else if (type=="interval") {
		#out <- .as.raster.interval(out, x)
	} else {
		out$range <- range
		#out <- .as.raster.continuous(out, x, type)
	}
	out$legend_type <- type
	if (!is.null(colNA)) {
		if (!is.na(colNA)) {
			out$colNA <- grDevices::rgb(t(grDevices::col2rgb(colNA)), alpha=alpha, maxColorValue=255)
			out$r[is.na(out$r)] <- out$colNA
		}
	}
	out
}


setMethod("plot", signature(x="SpatVector", y="character"), 
	function(x, y, col, type, mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, legend=TRUE, add=FALSE, plg=list(), pax=list(), ...)  {

		if (is.na(match(y, names(x)))) {
			stop(paste(y, "is not a name in x"))
		}
		if (missing(col)) {
			col <- topo.colors(100)
		}
	
		out <- .prep.vect.data(x, y, type=type, cols=col, mar=mar, draw=TRUE, plg=plg, pax=pax, legend=isTRUE(legend), axes=axes)
			
		out <- .plot.vect.map(x, out)
		invisible(out)
		
	}
)


setMethod("plot", signature(x="SpatVector", y="numeric"), 
	function(x, y, col, type, mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, legend=TRUE, add=FALSE, plg=list(), pax=list(), ...)  {
		if (y < 1) {
			plot(x, col=col, axes=axes, add=add, ...)
		} else if (y > ncol(x)) {
			stop(paste("x only has", ncol(x), " columns"))
		}
		plot(x, y=names(x)[y], cols=col, type=type, mar=mar, axes=axes, legend=legend, add=add, plg = plg, pax=pax, ...)
	}
)



setMethod("plot", signature(x="SpatVector", y="missing"), 
	function(x, y, col, axes=TRUE, add=FALSE, pax=list(), ...)  {
		if (missing(col)) {
			if (geomtype(x) == "points") {
				col <- "black"
			} else {
				col <- NULL
			}
		}
		col <- .getCols(size(x), col)
		out <- .prep.vect.data(x, y="", type="none", cols=col, mar=NULL, draw=TRUE, plg=list(), pax=pax, legend=FALSE, axes=axes)

		
		out <- .plot.vect.map(x, out)
		invisible(out)
	

		#
		#.vplot(x, NULL, col=col, axes=axes, add=add, ...)
	}
)


