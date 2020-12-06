# Author: Robert J. Hijmans
# Date :  June 2019
# Version 1.0
# License GPL v3


.plotLines <- function(x, out, lty=1, lwd=1, ...) {
	cols <- out$cols
	if (is.null(cols)) cols = rep("black", size(x))

	g <- geom(x)
	g <- split(g, g[,1])
	g <- lapply(g, function(x) split(x, x[,2]))
	#p <- sapply(g, function(x) lapply(x, function(y) lines(y[,3:4], ...))	
	n <- length(g)
	lty <- rep_len(lty, n)
	lwd <- rep_len(lwd, n)
	for (i in 1:n) {
		x <- g[[i]]
		for (j in 1:length(x)) {
			lines(x[[j]][,3:4], col=cols[i], lwd=lwd[i], lty=lty[i], ...)
		}
	}
	out$leg$lwd <- lwd
	out$leg$lty <- lty
	out
}

.plotPolygons <- function(x, out, lty=1, lwd=1, ...) {

	g <- geom(x)
	g <- split(g, g[,1])
	g <- lapply(g, function(y) split(y, y[,2]))
	n <- length(g)
	if (!is.null(out$leg$border)) {
		out$leg$border <- rep_len(out$leg$border, n)
	} else {
		out$leg$border <- NA
	}
	if (!is.null(out$leg$density)) {
		out$leg$density <- rep_len(out$leg$density, length(g))
		out$leg$angle <- rep_len(out$leg$angle, n)
	}
	out$leg$lty <- rep_len(lty, n)
	out$leg$lwd <- rep_len(lwd, n)

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
			if (!is.null(out$leg$density)) {
				graphics::polygon(a[,3:4], col=out$main_cols[i], density=out$leg$density[i], angle=out$leg$angle[i], border=NA, lwd=out$leg$lwd[i], lty=out$leg$lty[i], ...)
				graphics::polypath(a[,3:4], col=NA, rule="evenodd", border=out$border[i], lwd=out$leg$lwd[i], lty=out$leg$lty[i], ...)
			} else {
				graphics::polypath(a[,3:4], col=out$main_cols[i], rule = "evenodd", border=out$border[i], lwd=out$leg$lwd[i], lty=out$leg$lty[i], ...)
			}
		}
	}
	out
}


.vplot <- function(x, out, xlab="", ylab="", cex=1, pch=1, ...) {
	col <- out$main_cols
	if (out$leg$geomtype == "points") {
		#if (is.null(col)) col = "black"
		if (out$add) {
			points(x, col=col, cex=cex, pch=pch, ...)			
		} else {
			e <- as.vector(ext(x))
			plot(e[1:2], e[3:4], type="n", axes=FALSE, xlab=xlab, ylab=ylab, asp=out$asp)
			points(x, col=col, cex=cex, pch=pch, ...)			
		}
	} else {
		e <- matrix(as.vector(ext(x)), 2)
		if (out$leg$geomtype == "polygons") {
			out <- .plotPolygons(x, out, ...)
		} else {
			out <- .plotLines(x, out, ...)
		}
	}
	out
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


.vect.legend.none <- function(out) {
	if (out$leg$geomtype == "points") {
		out$main_cols <- .getCols(out$ngeom, out$cols)
	} else {
		out$cols <- .getCols(out$ngeom, out$cols)
	}
	out
}

.vect.legend.classes <- function(out) {
	ucols <- .getCols(length(out$uv), out$cols)
	uv <- sort(out$uv)
	i <- match(out$v, out$uv)
	out$cols <- ucols[i]
	out$main_cols <- out$cols

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
	
	if (!is.numeric(out$v)) {
		out$v <- as.integer(as.factor(out$v))
		z <- stats::na.omit(out$v)
		n <- length(z)
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

	brks <- seq(min(out$v, na.rm=TRUE), max(out$v, na.rm=TRUE), length.out = length(out$cols))
	grps <- cut(out$v, breaks = brks, include.lowest = TRUE)
	out$main_cols <- out$cols[grps]

	out
}


.vect.legend.interval <- function(out, ...) {

	nmx <- length(out$uv)
	if (nmx <= 1) {
		return(.vect.legend.classes(out, ...))
	}
	if (is.null(out$breaks)) {
		out$breaks <- 5
	} 
	out$breaks <- min(out$breaks, nmx)
	
	if (!is.numeric(out$v)) {
		out$v <- as.integer(as.factor(out$v))
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

	out$main_cols <- out$cols[out$vcut]
	out
}



.plot.vect.map <- function(x, out, xlab="", ylab="", type = "n", yaxs="i", xaxs="i", asp=out$asp, new=NA, density=NULL, angle=45, ...) {
	
	if ((!out$add) & (!out$legend_only)) {
		if (!any(is.na(out$mar))) { graphics::par(mar=out$mar) }
		plot(out$lim[1:2], out$lim[3:4], type="n", xlab=xlab, ylab=ylab, asp=asp, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
	}

	out$leg$density <- density
	out$leg$angle <- angle

	if (out$legend_type == "none") {
		out <- .vect.legend.none(out)
	} else if (out$legend_type == "classes") {
		out <- .vect.legend.classes(out)
	} else if (out$legend_type == "interval") {
		out <- .vect.legend.interval(out)
	} else {
		out <- .vect.legend.continuous(out)
		out$leg$density <- NULL
	}				
	if (!out$legend_only) {
		out <- .vplot(x, out, ...) 
	}

	if (out$axes) {
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


.prep.vect.data <- function(x, y, type, cols=NULL, mar=NULL, legend=TRUE, 
	legend.only=FALSE, levels=NULL, add=FALSE, range=NULL, new=NA, breaks=NULL, 
	xlim=NULL, ylim=NULL, colNA=NA, alpha=NULL, axes=TRUE, border="black",
	pax=list(), plg=list()) {

	out <- list()
	out$ngeom <- nrow(x)
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
	out$axs <- pax 
	out$leg <- plg
	out$leg$geomtype <- geomtype(x)
	out$asp <- 1
	out$lonlat <- isLonLat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(out$ext[3:4]) * pi)/180)
	}
	out$breaks <- breaks
	
	out$v <- unlist(x[, y, drop=TRUE], use.names=FALSE)
	out$uv <- sort(unique(out$v))

	if (missing(type)) {
		if (!is.numeric(out$uv) | length(out$uv) < 10) {
			type <- "classes"
		} else {
			type <- "continuous"
		}
	} else {
		type <- match.arg(type, c("continuous", "classes", "interval", "depends", "none"))
	}
	if (type=="none") {
		legend <- FALSE
		legend_only <- FALSE
	} else if (type=="classes") {
		out$levels <- levels
	} else if (type=="continuous") {
		out$range <- range
	}
	out$legend_type <- type

	if (missing(cols)) {
		if (type == "none") {
			if (out$leg$geomtype == "points") {
				cols <- "black"
			} else {
				cols <- NULL
			}
		} else {
			cols <- topo.colors(100)			
		}
	}
	if (!is.null(alpha)) {
		alpha <- clamp(alpha[1]*255, 0, 255)
		cols <- grDevices::rgb(t(grDevices::col2rgb(cols)), alpha=alpha, maxColorValue=255)
	} else {
		alpha <- 255
	}
	out$cols <- cols
	out$leg$border <- border
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)

	if (is.null(mar)) {
		if (out$legend_draw) {
			mar=c(3.1, 3.1, 2.1, 7.1)
		} else {
			mar=c(3.1, 3.1, 2.1, 2.1)
		}
	}
	out$mar <- mar

	if (!is.null(colNA)) {
		if (!is.na(colNA)) {
			out$colNA <- grDevices::rgb(t(grDevices::col2rgb(colNA)), alpha=alpha, maxColorValue=255)
			out$r[is.na(out$r)] <- out$colNA
		}
	}
	out
}


setMethod("plot", signature(x="SpatVector", y="character"), 
	function(x, y, col, type, mar=NULL, legend=TRUE, add=FALSE, axes=!add, main=y, plg=list(), pax=list(), nr, nc, ...) {

		y <- trimws(y)
		if (any(is.na(match(y, c("", names(x)))))) {
			i <- is.na(match(y, names(x)))
			stop(paste(paste(y[i], collapse=",")), "is not a name in x")
		}
		nrnc <- c(1,1)
		if (length(y) > 1) {
			nrnc <- .get_nrnc(nr, nc, length(y))
			old.par <- graphics::par(no.readonly =TRUE)
			on.exit(graphics::par(old.par))   
			graphics::par(mfrow=nrnc)
		}

		for (i in 1:length(y)) {
			if (length(y) > 1) {
				newrow <- (i %% nrnc[2]) == 1 
				lastrow <- i > (prod(nrnc) - nrnc[2])
				if (lastrow) {
					if (newrow) {
						pax$sides <- 1:2
					} else {
						pax$sides <- 1				
					}
				} else if (newrow) {
					pax$sides <- 2
				} else {
					pax$sides <- 0
				}
			}			
			if (y[i] == "") {
				out <- .prep.vect.data(x, y="", type="none", cols=col, mar=mar, plg=list(), pax=pax, legend=FALSE, add=add, axes=axes)
			} else {
				out <- .prep.vect.data(x, y[i], type=type, cols=col, mar=mar, plg=plg, pax=pax, legend=isTRUE(legend), add=add, axes=axes)
			}
			out <- .plot.vect.map(x, out, main=main[i], ...)
			invisible(out)		
		}
	}
)


setMethod("plot", signature(x="SpatVector", y="numeric"), 
	function(x, y, ...)  {
		y <- round(y)
		if (any(y > ncol(x))) {
			stop(paste("x only has", ncol(x), " columns"))
		}
		y[y<0] <- 0
		y <- c("", names(x))[y+1]
		plot(x, y, ...)
	}
)


setMethod("plot", signature(x="SpatVector", y="missing"), 
	function(x, y, ...)  {
		plot(x, "", ...)
	}
)

