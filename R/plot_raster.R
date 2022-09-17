
.as.raster.continuous <- function(out, x, type) {

	Z <- as.matrix(x, wide=TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA

# loss of precision
#	z <- stats::na.omit(round(as.vector(Z), 12))
	z <- stats::na.omit(as.vector(Z))
	n <- length(z)
	if (n == 0) {
		out$values = FALSE
		return(out)
	}

	uzi <- unique(z)

	if (type == "depends") {
		if (length(uzi) < 6) {
			return (.as.raster.classes(out, x))
		}
	} else if (length(uzi) == 1) {
		return (.as.raster.classes(out, x))
	}

	if (is.null(out$range)) {
		out$range <- range(z)
	} else {
		stopifnot(length(out$range) == 2)
		stopifnot(out$range[2] > out$range[1])
	}

	breaks <- .get_breaks(Z, length(out$cols), "eqint", out$range)
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

	Z <- as.matrix(x, wide=TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	if (all(is.na(Z))) {
		out$values = FALSE
		return(out)
	}

	fz <- as.factor(Z)
	levs <- as.integer(levels(fz))
	if (!is.null(out$levels)) {
		if (is.null(out$leg$legend)) {
			out$leg$legend <- as.character(out$levels)
		}
		levs <- out$levels
	} else {
		out$levels <- as.numeric(levs)
		if (is.null(out$leg$legend)) {
			out$leg$legend <- levs
		}
	}
	stopifnot(length(out$leg$legend) == length(out$levels))
	nlevs <- length(levs)

	if (NCOL(out$cols) == 2) {
		out$cols[,2] <- as.character(out$cols[,2])
		i <- match(Z, as.numeric(levs))
		Z[] <- out$cols[,2][i]
		i <- match(as.numeric(levs), out$cols[,1])
		out$leg$fill <- out$cols[i,2]
	} else {
		ncols <- length(out$cols)
		if (nlevs < ncols) {
			i <- round(seq(1, ncols, length.out = nlevs))
			cols <- out$cols[i]
		} else {
			cols <- rep_len(out$cols, nlevs)
		}
		out$leg$fill <- cols
		Z[] <- cols[as.numeric(fz)]
	}
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


.as.raster.factor <- function(out, x, ...) {

	z <- round(values(x))
	#z[z<1 | z>256] <- NA
	z[is.nan(z) | is.infinite(z)] <- NA
	if (all(is.na(z))) {
		out$values = FALSE
		return(out)
	}
	out$levels <- sort(stats::na.omit(unique(z)))
	ilevels <- match(out$levels, out$cats[[1]])
	# mismatch: ilevels <- na.omit(ilevels)
	if (out$all_levels) {
		out$leg$legend <- unique(na.omit(out$cats[, 2]))	
	} else {
		out$leg$legend <- unique(na.omit(out$cats[ilevels, 2]))
	}
	if (!is.null(out$coltab)) {
		if (out$all_levels) {
			mi <- match(out$cats[[1]], out$coltab[,1])
			mi[is.na(mi)] <- 1
			mc <- coltab[mi, ,drop=FALSE]
			out$leg$fill <- grDevices::rgb(mc[,2], mc[,3], mc[,4], mc[,5], maxColorValue=255)
		}
		out$levels <- out$levels[!is.na(ilevels)]
		m <- na.omit(match(out$cats[[1]][ilevels], out$coltab[,1]))
		out$coltab <- out$coltab[m, ,drop=FALSE]
		out$cols <- grDevices::rgb(out$coltab[,2], out$coltab[,3], out$coltab[,4], out$coltab[,5], maxColorValue=255)
		i <- match(z, out$coltab[,1])
		z <- out$cols[i]
	} else {
		levlab <- data.frame(id=out$levels, lab=out$cats[ilevels, 2], stringsAsFactors=FALSE)
		leglevs <- na.omit(unique(levlab[,2]))
		if (length(leglevs) == 0) {
			error("plot", "something is wrong with the categories")
		}
		nlevs <- length(leglevs)
		ncols <- length(out$cols)
		ncats <- nrow(out$cats)
		if (ncats < ncols) {
			i <- round(seq(1, ncols, length.out = ncats))
			out$cols <- out$cols[i]
		} else if (ncats > ncols) {
			out$cols <- rep_len(out$cols, ncats)
		}
		out$leg$fill <- out$cols
		out$cols <- out$cols[ilevels]
		dd <- data.frame(lab=leglevs, out$cols)
		m <- merge(levlab, dd)
		z <- m$out.cols[match(z, m$id)]
	}
	if (!out$all_levels) {
		out$leg$fill <- out$cols	
	}


	z <- matrix(z, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
	out$r <- as.raster(z)

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
	if (!is.null(out$legend_order)) {
		ord <- match(out$legend_order, out$leg$legend)
		out$leg$legend <- out$leg$legend[ord]
		out$leg$fill <- out$leg$fill[ord]
	} else if (isTRUE(out$legend_sort)) {
		ord <- order(out$leg$legend, decreasing=out$legend_sort_decreasing)
		out$leg$legend <- out$leg$legend[ord]
		out$leg$fill <- out$leg$fill[ord]
	}
	out
}


# to be merged with the vector variant.
.generic.interval <- function(out, Z) {
	if (is.null(out$breaks)) {
		out$breaks <- 5
	}
	if (length(out$breaks) == 1) {
		out$breaks <- .get_breaks(Z, out$breaks, out$breakby, out$range)
	}
	fz <- cut(Z, out$breaks, include.lowest=TRUE, right=FALSE)
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
	#out$cols <- cols
	out$leg$fill <- cols
	#out$leg$levels <- levels(fz)

	if (!is.null(out$leg$legend)) {
		stopifnot(length(out$leg$legend) == nlevs)
	} else {
		levs <- gsub("]", "", gsub(")", "", gsub("\\[", "", levs)))
		levs <- paste(levs, collapse=",")
		m <- matrix(as.numeric(unlist(strsplit(levs, ","))), ncol=2, byrow=TRUE)
		m <- apply(m, 1, function(i) paste(i, collapse=" - "))
		out$leg$legend <- m
	}
	out
}

.as.raster.interval <- function(out, x, ...) {

	out$legend_type <- "classes"

	if (NCOL(out$cols) == 3) {
		out$cols[,3] <- as.character(out$cols[,3])
		rcl <- cbind(as.matrix(out$cols[,1:2]), 1:nrow(out$cols))
		x <- classify(x, rcl, include.lowest=TRUE, others=NA)
		m <- apply(out$cols[,1:2], 1, function(i) paste(i, collapse=" - "))
		out$leg$legend <- m
		out$leg$fill <- out$cols[,3]
		Z <- as.matrix(x, wide=TRUE)
		Z[is.nan(Z) | is.infinite(Z)] <- NA
		Z[] <- out$leg$fill[Z]
	} else {
		Z <- as.matrix(x, wide=TRUE)
		Z[is.nan(Z) | is.infinite(Z)] <- NA
		out <- .generic.interval(out, Z)
		Z[] <- out$leg$fill[out$vcut]
	}
	if (is.null(out$leg$x)) { # && is.null(out$leg$ext)) {
		out$leg$x <- "top"
	}
	out$r <- as.raster(Z)
	out
}

# leg.shrink=c(0,0), leg.main=NULL, leg.main.cex = 1, leg.digits=NULL, leg.loc=NULL, leg.ext=NULL, leg.levels=NULL, leg.labels=NULL, leg.at=NULL,


.as.raster.colortable <- function(out, x, ...) {
	z <- round(values(x))
	#z[z<0 | z>255] <- NA
	z[is.nan(z) | is.infinite(z)] <- NA
	if (all(is.na(z))) {
		out$values = FALSE
		return(out)
	}
	out$cols <- grDevices::rgb(out$coltab[,2], out$coltab[,3], out$coltab[,4], out$coltab[,5], maxColorValue=255)
	i <- match(z, out$coltab[,1])
	z <- out$cols[i]
	z <- matrix(z, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
	out$r <- as.raster(z)
	out
}


# legend
#	border="black", box.lwd = graphics::par("lwd"), box.lty = graphics::par("lty"),
#	box.col = graphics::par("fg"), bty = "o", bg = graphics::par("bg"), xjust = 0,
#	yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), text.width = NULL,
#	text.col = graphics::par("col"), text.font = NULL, ncol = 1, horiz = FALSE, title = NULL,
 #   inset = 0, title.col = text.col, title.adj = 0.5,

.plotit <- function(x, xlab="", ylab="", type = "n", yaxs="i", xaxs="i", asp=x$asp, axes=TRUE, new=NA,
	main="", line=0.5, cex.main=0.8, font.main=graphics::par()$font.main,
	col.main = graphics::par()$col.main, reset=FALSE, ...) {

#	if (x$add) axes = FALSE
	if ((!x$add) & (!x$legend_only)) {

		old.mar <- graphics::par()$mar
		if (!any(is.na(x$mar))) { graphics::par(mar=x$mar) }
		if (reset) on.exit(graphics::par(mar=old.mar))

		plot(x$lim[1:2], x$lim[3:4], type=type, xlab=xlab, ylab=ylab, asp=asp, xaxs=xaxs, yaxs=yaxs, axes=!x$values, ...)

		graphics::title(main, line=line, cex.main=cex.main, font.main=font.main, col.main=col.main)
	}
	if (!x$values) {
		return(x)
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
			leg <- do.call(.plot.class.legend, x$leg)
		}
	}
	x
}



.prep.plot.data <- function(x, type, cols, mar=NULL, draw=FALSE,
  interpolate=FALSE, legend=TRUE, legend.only=FALSE, pax=list(), plg=list(),
  levels=NULL, add=FALSE, range=NULL, new=NA, breaks=NULL, breakby="eqint",
  coltab=NULL, cats=NULL, xlim=NULL, ylim=NULL, ext=NULL, colNA=NA, alpha=NULL, reset=FALSE,
  sort=TRUE, decreasing=FALSE, grid=FALSE, las=0, all_levels=FALSE, ...) {

#mar=c(5.1, 4.1, 4.1, 7.1); legend=TRUE; axes=TRUE; pal=list(); pax=list(); maxcell=50000; draw=FALSE; interpolate=FALSE; legend=TRUE; legend.only=FALSE; pax=list(); pal=list(); levels=NULL; add=FALSE; range=NULL; new=NA; breaks=NULL; coltab=NULL; facts=NULL; xlim=NULL; ylim=NULL;

	out <- list()

	out$lim <- out$ext <- as.vector(ext(x))
	if (!is.null(ext)) {
		ext <- ext(ext)
		x <- crop(x, ext)
		out$ext <- as.vector(ext(x))
		out$lim <- ext
	}
	if (!(is.null(xlim) & is.null(ylim))) {
		e <- as.vector(ext(x))
		if (!is.null(xlim)) e[1:2] <- sort(xlim)
		if (!is.null(ylim)) e[3:4] <- sort(ylim)
		x <- crop(x, ext(e))
		out$ext <- as.vector(ext(x))
		out$lim <- e
	}

	out$add <- isTRUE(add)
	out$axs <- as.list(pax)
	if (is.null(out$axs$las)) out$axs$las <- las
	out$draw_grid <- isTRUE(grid)

	out$leg <- as.list(plg)
	out$all_levels <- isTRUE(all_levels)

	out$asp <- 1
	out$lonlat <- is.lonlat(x, perhaps=TRUE, warn=FALSE)
	if (out$lonlat) {
		out$asp <- 1/cos((mean(out$ext[3:4]) * pi)/180)
	}

	if (!is.null(alpha)) {
		if (!inherits(alpha, "SpatRaster")) {
			cols <- grDevices::rgb(t(grDevices::col2rgb(cols)), alpha=alpha[1]*255, maxColorValue=255)
		}
	} else {
		alpha <- 255
	}

	out$cols <- cols
	out$coltab <- coltab
	out$cats <- cats
	out$breaks <- breaks
	out$breakby <- breakby
	out$interpolate <- FALSE
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)
	if (!is.logical(sort)) {
		out$legend_order <- sort
		out$legend_sort <- FALSE
	} else {
		out$legend_sort <- isTRUE(sort)
	}
	out$legend_sort_decreasing <- isTRUE(decreasing)

	if (is.null(mar)) {
		if (out$legend_draw) {
			out$mar <- c(3.1, 3.1, 2.1, 7.1)
		} else {
			out$mar <- c(3.1, 3.1, 2.1, 2.1)
		}
	} else {
		out$mar <- rep_len(mar, 4)
	}

	if (!hasValues(x)) {
		out$values = FALSE
		warn("plot", "SpatRaster has no cell values")
	} else {
		out$values <- TRUE

		if (type=="factor") {
			out <- .as.raster.factor(out, x)
		} else if (type=="colortable") {
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
			if (!is.na(colNA) && out$values) {
				out$colNA <- grDevices::rgb(t(grDevices::col2rgb(colNA)), alpha=alpha, maxColorValue=255)
				out$r[is.na(out$r)] <- out$colNA
			}
		}
	}

	if (draw) {
		if (inherits(alpha, "SpatRaster")) {
			alpha <- clamp(as.vector(alpha[[1]])*255, 0, 255)
			out$r <- matrix(grDevices::rgb(t(grDevices::col2rgb(out$r)), alpha=alpha, maxColorValue=255),
			nrow=nrow(out$r), byrow=TRUE)
		}

		out <- .plotit(out, new=new, reset=reset, ...)
	}
	invisible(out)
}


setMethod("plot", signature(x="SpatRaster", y="numeric"),
	function(x, y=1, col, type, mar=NULL, legend=TRUE, axes=TRUE, plg=list(), pax=list(), maxcell=500000, smooth=FALSE, range=NULL, levels=NULL, all_levels=FALSE, breaks=NULL, breakby="eqint", fun=NULL, colNA=NULL, alpha=NULL, sort=FALSE, decreasing=FALSE, grid=FALSE, ext=NULL, reset=FALSE, ...) {

		y <- round(y)
		stopifnot((min(y) > 0) & (max(y) <= nlyr(x)))

		if (length(y) > 1) {
			x <- x[[y]]
			if (inherits(alpha, "SpatRaster")) {
				if (nlyr(alpha) > 1) {
					alpha <- alpha[[y]]
				}
			}
			plot(x, col=col, type=type, mar=mar, legend=legend, axes=axes, plg=plg, pax=pax, maxcell=maxcell/(length(x)/2), smooth=smooth, range=range, levels=levels, all_levels=all_levels, breaks=breaks, breakby=breakby, fun=fun, colNA=colNA, alpha=alpha, grid=grid, sort=sort, decreasing=decreasing, ext=ext, ...)
			return(invisible())
		}

		if (inherits(alpha, "SpatRaster")) {
			if (!compareGeom(x, alpha, crs=FALSE, ext=FALSE, rowcol=TRUE)) {
				error("plot", "geometry of alpha does not match x")
			}
		}

		x <- x[[y]]
		if (ncell(x) > 1.1 * maxcell) {
			if (inherits(alpha, "SpatRaster")) {
				if (nlyr(alpha) > 1) {
					alpha <- alpha[[y]]
				}
				alpha <- spatSample(alpha, maxcell, ext=ext, method="regular", as.raster=TRUE)
			}
			x <- spatSample(x, maxcell, ext=ext, method="regular", as.raster=TRUE)
		}
		
		if (is.character(legend)) {
			plg$x <- legend
			legend <- TRUE
		}

		if (missing(col)) {
			col <- rev(grDevices::terrain.colors(255))
		} else if (inherits(col, "data.frame")) {
			if (ncol(col) == 2) {
				type <- "classes"
			} else if (ncol(col) == 3) {
				type <- "interval"
			} else {
				error("plot", "number of columns of a col data.frame should be 2 or 3")
			}
			breaks <- NULL
		}
		coltab <- NULL
		cats  <- NULL
		if (!is.null(breaks)) {
			if (missing(type)) {
				type <- "interval"
			} else {
				range <- range(breaks)
			}
		} else {
			if (missing(type)) {
				if (has.colors(x)) {
					coltab <- coltab(x)[[1]]
					if (is.factor(x)) {
						#act <- activeCat(x)
						cats <- levels(x)[[1]] # cats(x)[[1]][, c(1, act+1)]
						type <- "factor"
					} else {
						type <- "colortable"
						legend <- FALSE
					}
				} else if (is.factor(x)) {
					type <- "factor"
					#act <- activeCat(x)
					cats <- levels(x)[[1]] #cats(x)[[1]][, c(1, act+1)]
				} else if (is.bool(x)) {
					type <- "factor"
					levels(x) <- data.frame(id=0:1, value=c("False", "True"))
					cats <- cats(x)[[1]]
				} else {
					type <- "depends"
				}
			} else {
				type <- match.arg(type, c("continuous", "classes", "interval"))
			}
		}

		x <- .prep.plot.data(x, type=type, cols=col, mar=mar, draw=TRUE, plg=plg, pax=pax, legend=isTRUE(legend), axes=isTRUE(axes), coltab=coltab, cats=cats, interpolate=smooth, levels=levels, range=range, colNA=colNA, alpha=alpha, reset=reset, grid=grid, sort=sort, decreasing=decreasing, ext=ext, all_levels=all_levels, breaks=breaks, breakby=breakby, ...)

		if (!is.null(fun)) {
			if (!is.null(formals(fun))) {
				fun(y)
			} else {
				fun()
			}
		}
		invisible(x)
	}
)



setMethod("plot", signature(x="SpatRaster", y="missing"),
	function(x, y, maxcell=500000, main, mar=NULL, nc, nr, maxnl=16, legend, ...)  {

		
		if (has.RGB(x)) {
			i <- x@ptr$getRGB() + 1
			if (missing(main)) main = ""
			if (is.null(mar)) mar = 0
			if (missing(legend)) {
				legend <- FALSE
			} else {
				legend <- legend
			}
			plotRGB(x, i[1], i[2], i[3], maxcell=maxcell, mar=mar, main=main, ...)
			return(invisible())
		}

		nl <- max(1, min(nlyr(x), maxnl))

		if (missing(legend)) legend <- TRUE
		
		if (nl==1) {
			if (missing(main)) main = ""
			out <- plot(x, 1, maxcell=maxcell, main=main[1], mar=mar, legend=legend, ...)
			return(invisible(out))
		}

		nrnc <- .get_nrnc(nr, nc, nl)
		old.par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old.par))
		if (is.null(mar)) {
			mar=c(1.5, 1.5, 1.5, 3)
		}
		graphics::par(mfrow=nrnc)
		maxcell=maxcell/(nl/2)

		if (missing("main")) {
			tm <- time(x)
			if (!any(is.na(tm))) {
				main <- as.character(time(x))
			} else {
				main <- names(x)
			}
		} else {
			main <- rep_len(main, nl)
		}
		#if (onelegend) { legend <- FALSE }
		legend <- rep_len(legend, length.out=nl)
		for (i in 1:nl) {
			plot(x, i, main=main[i], mar=mar, legend=legend[i], ...)
		}
	}
)


setMethod("plot", signature(x="SpatRaster", y="character"),
	function(x, y, ...) {
		y <- match(y, names(x))
		if (any(is.na(y))) {
			error("plot", "y does not match the names in x")
		}
		plot(x, y, ...)
	}
)

