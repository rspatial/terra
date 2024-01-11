
.as.raster.rgb <- function(out, x) {

	if (is.null(out$rgb$scale)) {
		scale <- 255
		if ( all(hasMinMax(x)) ) {
			rng <- minmax(x)[, 1:3]
			scale <- max(max(rng[2]), 255)
		}
	} else {
		scale <- out$rgb$scale
	}
	
	if (!is.null(out$rgb$zlim)) {
		if (length(out$rgb$zlim) == 2) {
			out$rgb$zlim <- sort(out$rgb$zlim)
			if (isTRUE(out$rgb$zcol)) {
				x <- clamp(x, out$rgb$zlim[1], out$rgb$zlim[2], values=TRUE)
			} else { #if (is.na(zlimcol)) {
				x <- clamp(x, out$rgb$zlim[1], out$rgb$zlim[2], values=FALSE)
			}
		} else if (NROW(out$rgb$zlim) == 3 & NCOL(out$rgb$zlim) == 2) {
			for (i in 1:3) {
				zmin <- min(out$rgb$zlim[i,])
				zmax <- max(out$rgb$zlim[i,])
				if (isTRUE(out$rgb$zcol)) {
					x[[i]] <- clamp(x[[i]], zmin, zmax, values=TRUE)
				} else { #if (is.na(zlimcol)) {
					x[[i]] <- clamp(x[[i]], zmin, zmax, values=FALSE)
				}
			}
		} else {
			error('zlim should be a vector of two numbers or a 3x2 matrix (one row for each color)')
		}
	}

	if (!is.null(out$rgb$stretch)) {
		if (out$rgb$stretch == "lin") {
			if ((!is.null(out$rgb$zlim)) && (length(out$rgb$zlim) == 2)) {
				x <- stretch(x, smin=out$rgb$zlim[1], smax=out$rgb$zlim[2])
			} else {
				x <- stretch(x, minq=0.02, maxq=0.98)
			}
		} else {
			x <- stretch(x, histeq=TRUE, scale=255)
		}
		scale <- 255
	}

	RGB <- values(x)
	RGB <- stats::na.omit(RGB)
	naind <- as.vector( attr(RGB, "na.action") )

	if (ncol(RGB) == 4){
		alpha <- RGB[,4] * 255
		RGB <- RGB[,-4]
	}

	if (!is.null(naind)) {
		bg <- grDevices::col2rgb(out$rgb$colNA)
		if (is.null(out$rgb$bgalpha)) out$rgb$bgalpha <- 255
		bg <- grDevices::rgb(bg[1], bg[2], bg[3], alpha=out$rgb$bgalpha, maxColorValue=255)
		z <- rep( bg, times=ncell(x))
		z[-naind] <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=out$alpha, maxColorValue=scale)
	} else {
		z <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=out$alpha, maxColorValue=scale)
	}
	
	out$r <- matrix(z, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
	out
}


.as.raster.continuous <- function(out, x, type) {

	Z <- as.matrix(x, wide=TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA

# loss of precision
#	z <- stats::na.omit(round(as.vector(Z), 12))
	z <- stats::na.omit(as.vector(Z))
	n <- length(z)
	if (n == 0) {
		#out$values = FALSE
		out$range <- c(NA, NA)
		out$legend_draw <- FALSE
		return(out)
	}

	uzi <- round(unique(z), 12)

	if (type == "depends") {
		if (length(uzi) < 9) {
			return (.as.raster.classes(out, x))
		}
	} else if ((length(uzi) == 1) && is.null(out$range)) {
		return (.as.raster.classes(out, x))
	}

	if (is.null(out$range)) {
		out$range <- range(z)
#		out$fill_range <- FALSE
	} else {
		stopifnot(length(out$range) == 2)
		stopifnot(out$range[2] > out$range[1])
	}

	breaks <- .get_breaks(z, length(out$cols), "eqint", out$range)
#	if (isTRUE(out$fill_range)) {
#		zrng <- range(z)
#		breaks[1] <- zrng[1]
#		breaks[length(breaks)] <- zrng[2]
#		out$frange <- zrng
#	} else {
#		out$frange <- out$range	
#	}
	Z[] <- out$cols[as.integer(cut(Z, breaks, include.lowest=TRUE, right=FALSE))]

	out$r <- as.raster(Z)

	out$legend_type <- "continuous"

	if (is.null(out$levels)) {
		out$levels <- 5
	}
	if (is.null(out$leg$digits)) {
		dif <- diff(out$range)
		if ((dif == 0) || (length(dif) ==0)) {
			out$leg$digits <- 0;
		} else {
			out$leg$digits <- max(0, -floor(log10(dif/10)))
		}
	}
	if (is.null(out$leg$x)) out$leg$x <- "right"
	out
}


prettyNumbs <- function(x, digits) {
	x <- formatC(x, digits=digits, format = "f", flag="#")
	x <- substr(x, 1, digits+1)
	gsub("\\.$", "", x)
}

.as.raster.classes <- function(out, x, ...) {

	Z <- as.matrix(x, wide=TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	if (all(is.na(Z))) {
		#out$values = FALSE
		out$range <- c(NA, NA)
		out$legend_draw <- FALSE
		return(out)
	}

	fz <- as.factor(Z)
	if (!is.null(out$levels)) {
		if (is.null(out$leg$legend)) {
			out$leg$legend <- as.character(out$levels)
		}
		levs <- out$levels
	} else {
		levs <- as.numeric(levels(fz))
		digits <- out$leg$digits
		if (is.null(digits)) {
			if (length(levs) > 1) {
				d <- ceiling(1 / min(diff(sort(levs))))
				decimals <- round(log10(d) + 1)
			} else {
				txt <- format(levs, scientific = FALSE, digits=18)
				txt <- unlist(strsplit(txt, "\\."))
				if (nchar(txt[1]) > 3) decimals <- 0
				else if (nchar(txt[1]) > 2) decimals <- 1
				else if (nchar(txt[1]) > 1) decimals <- 2
				else if (txt[1] != "0") decimals <- 3
				else if (length(txt) > 1) {
					txt <- unlist(strsplit(txt[2], ""))
					i <- which(txt != "0")[1]
					if (length(i) > 0) decimals <- i+2
					else decimals <- 9;
				} else {
					decimals <- 0
				}
			}
			levs <- round(levs, decimals)
		}
		out$levels <- levs
		if (is.null(out$leg$legend)) {
			if (!is.null(out$leg$digits)) {
				out$leg$legend <- prettyNumbs(levs, digits)
			} else {
				out$leg$legend <- levs
			}
		}
	}
	out$leg$digits <- NULL
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
		if (nlevs == 1) {
			cols <- out$cols[length(out$cols)]
		} else if (nlevs < ncols) {
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
			out$leg$x = "default"
			out$leg$y = NULL
		} else {
			if (length(out$leg$ext) == 4) {
				out$leg$x = out$leg$ext[1]
				out$leg$y = out$leg$ext[4]
			} else {
				out$leg$x = "default"
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
		#out$values = FALSE
		out$range <- c(NA, NA)
		out$legend_draw <- FALSE
		return(out)
	}
	out$levels <- sort(stats::na.omit(unique(z)))
	if (out$all_levels) {
		ilevels <- 1:nrow(out$cats)
		out$levels <- out$cats[,1]
	} else {
		ilevels <- match(out$levels, out$cats[[1]])
		if (any(is.na(ilevels))) {
			warn("plot", "unknown categories in raster values")
		}
	}

	if (!is.null(out$coltab)) {
		# avoid multiple colors for the same category
		ilevs <- stats::na.omit(ilevels)
		ulevs <- unique(out$cats[ilevs,2])
		if (length(ulevs) < length(ilevs)) {
			z <- out$cats[match(z, out$cats[,1]),2]
			i <- match(ulevs, out$cats[ilevels,2])	
			j <- out$cats[ilevs[i],]
			z <- j[match(z, j[,2]), 1]
			out$cats <- j			
			out$levels <- sort(stats::na.omit(unique(z)))
			#out$coltabt = out$coltab[match(j[,1], out$coltab[,1]), ]
			if (out$all_levels) {
				ilevels <- 1:nrow(out$cats)
				out$levels <- out$cats[,1]
			} else {
				ilevels <- match(out$levels, out$cats[[1]])
			}
		}

		if (out$all_levels) {
			mi <- match(out$cats[[1]], out$coltab[,1])
			mi[is.na(mi)] <- 1
			mc <- out$coltab[mi, ,drop=FALSE]
			out$leg$fill <- grDevices::rgb(mc[,2], mc[,3], mc[,4], mc[,5], maxColorValue=255)
			if (is.null(out$leg$legend)) out$leg$legend <- stats::na.omit(out$cats[, 2])
		} else {	
			out$levels <- out$levels[!is.na(ilevels)]
			m <- stats::na.omit(match(out$cats[[1]][ilevels], out$coltab[,1]))
			if (is.null(out$leg$legend)) out$leg$legend <- stats::na.omit(out$cats[ilevels, 2])
			out$coltab <- out$coltab[m, ,drop=FALSE]
		}

		out$cols <- grDevices::rgb(out$coltab[,2], out$coltab[,3], out$coltab[,4], out$coltab[,5], maxColorValue=255)
		i <- match(z, out$coltab[,1])
		z <- out$cols[i]
	} else {
		if (is.null(out$leg$legend)) out$leg$legend <- unique(stats::na.omit(out$cats[ilevels, 2]))
		levlab <- data.frame(id=out$levels, lab=out$cats[ilevels, 2], stringsAsFactors=FALSE)
		leglevs <- stats::na.omit(unique(levlab[,2]))
		if (length(leglevs) == 0) {
			error("plot", "something is wrong with the categories")
		}
		nlevs <- length(leglevs)
		ncols <- length(out$cols)
		#ncats <- nrow(out$cats)
		if (nlevs < ncols) {
			i <- round(seq(1, ncols, length.out = nlevs))
			out$cols <- out$cols[i]
		} else if (nlevs > ncols) {
			out$cols <- rep_len(out$cols, nlevs)
		}
		out$leg$fill <- out$cols
		#out$cols <- out$cols[ilevels] 
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
			out$leg$x = "default"
			out$leg$y = NULL
		} else {
			if (length(out$leg$ext) == 4) {
				out$leg$x = out$leg$ext[1]
				out$leg$y = out$leg$ext[4]
			} else {
				out$leg$x = "default"
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

	if (!is.null(out$leg$digits)) {
#		out$leg$legend <- substr(formatC(levs, digits=digits, format = "f", flag="#"), 1, digits+1)
		fz <- cut(Z, out$breaks, include.lowest=TRUE, right=FALSE, dig.lab=out$leg$digits)
	} else {
		fz <- cut(Z, out$breaks, include.lowest=TRUE, right=FALSE)
	}


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
		if (!is.null(out$leg$digits)) {
			m <- prettyNumbs(m, out$leg$digits)
		}
		m <- apply(m, 1, function(i) paste(i, collapse=" - "))
		out$leg$legend <- m
	}
	out$leg$digits <- NULL
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
		out$leg$x <- "default"
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

.plotit <- function(x) {

	if (is.null(x$r)) {
		x$values = FALSE
	}

	if (x$add) {
		reset.clip()
	} else if (!x$legend_only) {
		old.mar <- graphics::par()$mar
		if (!any(is.na(x$mar))) { graphics::par(mar=x$mar) }
		if (x$reset) on.exit(graphics::par(mar=old.mar))
		arglist <- c(list(x=x$lim[1:2], y=x$lim[3:4], type="n", xlab="", ylab="", asp=x$asp, xaxs=x$xaxs, yaxs=x$yaxs, axes=FALSE), x$dots)
		do.call(plot, arglist)
		if (!is.null(x$background)) {
			graphics::rect(x$lim[1], x$lim[3], x$lim[2], x$lim[4], col=x$background)			
		}
	}
	if (!x$values) {
		if (!x$add) try(set.clip(x$lim, x$lonlat))
		return(x)
	}
	if (!x$legend_only) {
		graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4],
			angle = 0, interpolate = x$interpolate)

		if (x$axes) x <- .plot.axes(x)
	}
	
	if (x$legend_draw) {
		if (x$legend_type == "continuous") {
			x <- do.call(.plot.cont.legend, list(x=x))
#		} else if (x$legend_type == "classes") {
		} else {
			if (x$add) {
				if (x$clip) {
					x$leg$plotlim <- unlist(get.clip()[1:4])
				} else {
					x$leg$plotlim <- graphics::par("usr")
				}			
				if (is.null(x$leg$plotlim)) {
					x$leg$plotlim <- x$lim			
				}				
			} else {
				if (x$clip) {
					x$leg$plotlim <- x$lim
				} else {
					x$leg$plotlim <- graphics::par("usr")
				}				
			}
			x$leg$used <- do.call(.plot.class.legend, x$leg)
		}
	}
	if (isTRUE(x$box)) { 
		if (x$clip) {
			lines(ext(x$lim))	
		} else {
			lines(ext(graphics::par("usr")))		
		}
	}
	
	if ((x$main != "") && (!x$add) && (!x$legend_only)) {
		pos <- 3
		if (is.null(x$loc.main)) {
			if (isTRUE(x$clip)) {
				x$loc.main <- c(x$lim[1] + diff(x$lim[1:2]) / 2, x$lim[4])
			} else {
				usr <- graphics::par("usr")			
				x$loc.main <- c(usr[1] + diff(usr[1:2]) / 2, usr[4])			
			}
		} else if (inherits(x$loc.main, "character")) {
			xyp <- .txt.loc(x)
			x$loc.main <- xyp[1:2]
			pos <- xyp[3]
		}
		if (isTRUE(x$halo.main)) {
			.halo(x$loc.main[1], x$loc.main[2], x$main, pos=pos, offset=x$line.main, cex=x$cex.main, 
				font=x$font.main, col=x$col.main, xpd=TRUE)
		} else {
			text(x$loc.main[1], x$loc.main[2], x$main, pos=pos, offset=x$line.main, cex=x$cex.main, 
				font=x$font.main, col=x$col.main, xpd=TRUE)
		}
	}
	if (!x$add) {
		try(set.clip(x$lim, x$lonlat))
	}
	invisible(x)
}



.prep.plot.data <- function(x, type, cols, mar=NULL, draw=FALSE,
  interpolate=FALSE, legend=TRUE, legend.only=FALSE, pax=list(), plg=list(),
  levels=NULL, add=FALSE, range=NULL, breaks=NULL, breakby="eqint",
  coltab=NULL, cats=NULL, xlim=NULL, ylim=NULL, ext=NULL, colNA=NA, alpha=NULL, reset=FALSE,
  sort=TRUE, decreasing=FALSE, grid=FALSE, las=0, all_levels=FALSE, decimals=NULL, background=NULL,
  xlab="", ylab="", cex.lab=0.8, line.lab=1.5, asp=NULL, yaxs="i", xaxs="i", main="", cex.main=1.2, 
  line.main=0.5, font.main=graphics::par()$font.main, col.main = graphics::par()$col.main, loc.main=NULL, 
  halo=FALSE, axes=TRUE, box=TRUE, cex=1, maxcell=500000, buffer=FALSE, clip=TRUE, 
  # for rgb 
  stretch=NULL, scale=NULL, bgalpha=NULL, zlim=NULL, zcol=NULL,  ...) {
#cex is catch and kill

	out <- list()
	e <- out$lim <- out$ext <- as.vector(ext(x))
	hadWin <- hasWin <- FALSE
	if (add && is.null(ext)) {
		ext <- unlist(get.clip())[1:4]
	}
	
	if ((!is.null(ext)) || (!is.null(xlim)) || (!is.null(ylim))) {
		if (!is.null(ext)) {
			ext <- ext(ext)
			#e <- as.vector(align( intersect(ext, ext(x)), x))
			e <- as.vector(align(ext, x))
			out$lim <- out$ext <- e
		} 
		if (!is.null(xlim)) {
			stopifnot(length(xlim) == 2)
			e[1:2] <- sort(xlim)
		}
		if (!is.null(ylim)) {
			stopifnot(length(ylim) == 2)
			e[3:4] <- sort(ylim)
		}
		out$lim <- e
		
		hasWin <- TRUE
		hadWin <- window(x)
		oldWin <- ext(x)
		w <- intersect(ext(x), ext(e))		
		window(x) <- out$ext <- w
	} 
	if (ncell(x) > 1.1 * maxcell) {
		if (inherits(alpha, "SpatRaster")) {
			if (nlyr(alpha) > 1) {
				alpha <- alpha[[1]]
			}
			alpha <- spatSample(alpha, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
		}
		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
		out$lim <- out$ext <- as.vector(ext(x))
	}
	
	if (buffer) {
		dx <- diff(out$lim[1:2]) / 50
		dy <- diff(out$lim[3:4]) / 50
		out$lim[1:2] <- out$lim[1:2] + c(-dx, dx)
		out$lim[3:4] <- out$lim[3:4] + c(-dy, dy)
	}

	out$add <- isTRUE(add)
	out$axs <- as.list(pax)
	if (is.null(out$axs$las)) out$axs$las <- las
	if (is.null(out$axs$cex.lab)) out$axs$cex.lab <- cex.lab
	if (is.null(out$axs$line.lab)) out$axs$line.lab <- line.lab

	out$draw_grid <- isTRUE(grid)
	out$leg$digits <- decimals

	out$leg <- as.list(plg)
	out$all_levels <- isTRUE(all_levels)

	if (is.null(asp)) {
		out$lonlat <- is.lonlat(x, perhaps=TRUE, warn=FALSE)
		if (out$lonlat) {
			out$asp <- 1/cos((mean(out$ext[3:4]) * pi)/180)
		} else {
			out$asp <- 1
		}
	} else {
		out$asp <- asp
		out$lonlat <- FALSE
	}
	if (!is.null(alpha)) {
		if (!inherits(alpha, "SpatRaster")) {
			alpha <- alpha[1] * 255
			cols <- grDevices::rgb(t(grDevices::col2rgb(cols)), alpha=alpha, maxColorValue=255)
		} 
	} else {
		alpha <- 255
	}
	out$rgb$stretch <- stretch
	out$rgb$scale <- scale
	out$rgb$bgalpha <- bgalpha
	out$rgb$zlim <- zlim
	out$rgb$zcol <- isTRUE(zcol)
	if (is.null(colNA) || is.na(colNA)) {
		out$rgb$colNA <- "white"
	} else {
		out$rgb$colNA = colNA
	}
	out$clip <- isTRUE(clip)
	out$dots  <- list(...)
	out$reset <- reset
	out$main  <- main
	out$halo.main <- halo
	out$loc.main  <- loc.main
	out$cex.main  <- cex.main
	out$font.main <- font.main
	out$col.main  <- col.main
	out$line.main <- line.main
	out$axes <- axes
	out$xaxs <- xaxs
	out$yaxs <- yaxs
	out$xlab <- xlab
	out$ylab <- ylab 
	out$cols <- cols
	out$coltab <- coltab
	out$cats <- cats
	out$breaks <- breaks
	out$breakby <- breakby
	out$interpolate <- FALSE
	out$background <- background
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)
	if (!is.logical(sort)) {
		out$legend_order <- sort
		out$legend_sort <- FALSE
	} else {
		out$legend_sort <- isTRUE(sort)
	}
	out$legend_sort_decreasing <- isTRUE(decreasing)
	out$box <- isTRUE(box)
	
	if (!is.null(out$leg$loc)) {
		out$leg$x <- out$leg$loc
		out$leg$loc <- NULL
	}


	if (!hasValues(x)) {
		out$values <- FALSE
		out$legend_draw <- FALSE
		#warn("plot", "SpatRaster has no cell values")
	} else {
		out$values <- TRUE

		if (type=="factor") {
			out <- .as.raster.factor(out, x)
		} else if (type=="rgb") {
			out <- .as.raster.rgb(out, x)
			out$interpolate <- isTRUE(interpolate)
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
			#out$fill_range <- fill_range
			out <- .as.raster.continuous(out, x, type)
		}

		if (is.null(mar)) {
			out$mar <- c(2, 2, 2, 2)
			if (out$legend_draw) {
				if (is.null(out$leg$ext)) {
					if (is.null(out$leg$x)) {
						out$leg$x <- "default"
						out$mar <- c(2, 2, 2, 4)
					} else if (out$legend_type == "continuous") {
						if (out$leg$x == "top") {
							out$mar <- c(2, 2, 4, 2)
						} else if (out$leg$x == "bottom") {
							out$mar <- c(4, 2, 2, 2)
						} else if (out$leg$x == "left") {
							out$mar <- c(2, 5, 2, 1)
						} else {
							out$mar <- c(2, 2, 2, 4)
						}
					} else if (out$leg$x == "default") {
						out$mar <- c(2, 2, 2, 4)					
					}
				} 
			}
		} else {
			out$mar <- rep_len(mar, 4)
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
		out <- .plotit(out)
	}
	
	if (any(hasWin)) {
		window(x) <- NULL
		if (any(hadWin)) {
			window(x) <- oldWin
		}
	}
	invisible(out)
}


setMethod("plot", signature(x="SpatRaster", y="numeric"),
	function(x, y=1, col, type=NULL, mar=NULL, legend=TRUE, axes=!add, plg=list(), pax=list(), maxcell=500000, smooth=FALSE, range=NULL, levels=NULL, all_levels=FALSE, breaks=NULL, breakby="eqint", fun=NULL, colNA=NULL, alpha=NULL, sort=FALSE, decreasing=FALSE, grid=FALSE, ext=NULL, reset=FALSE, add=FALSE, buffer=FALSE, background=NULL, box=axes, clip=TRUE, ...) {

		old.mar <- graphics::par()$mar
		on.exit(graphics::par(mar=old.mar))

		
		y <- round(y)
		hasRGB <- FALSE		

		if (has.RGB(x) && ((is.null(type) && (y[1] < 0)))) {
			type <- "rgb"
			legend <- FALSE
			if (is.null(mar)) {
				mar <- 0
				axes <- FALSE
			}
			hasRGB <- TRUE
			y <- RGB(x)
		}
		nlx <- nlyr(x)
		stopifnot((min(y) > 0) & (max(y) <= nlx))

		if ((!hasRGB) && (length(y) > 1)) {

			x <- x[[y]]
			if (inherits(alpha, "SpatRaster")) {
				if (nlyr(alpha) > 1) {
					alpha <- alpha[[y]]
				}
			}
			plot(x, col=col, type=type, mar=mar, legend=legend, axes=axes, plg=plg, pax=pax, maxcell=2*maxcell/length(y), smooth=smooth, range=range, levels=levels, all_levels=all_levels, breaks=breaks, breakby=breakby, fun=fun, colNA=colNA, alpha=alpha, grid=grid, sort=sort, decreasing=decreasing, ext=ext, reset=reset, add=add, buffer=buffer, background=background, box=box, clip=clip, ...)
			return(invisible())
		} else {
			x <- x[[y]]
		}

		if (inherits(alpha, "SpatRaster")) {
			if (!compareGeom(x, alpha, crs=FALSE, ext=FALSE, rowcol=TRUE)) {
				error("plot", "geometry of alpha does not match x")
			}
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
			if (is.null(type)) {
				type <- "interval"
			} else {
				range <- range(breaks)
			}
		} else {
			if (is.null(type)) {
				if (has.colors(x)) {
					coltab <- coltab(x)[[1]]
					if (is.factor(x)) {
						if (activeCat(x) >= 0) {
							cats <- levels(x)[[1]] 
							type <- "factor"
						} else {
							type <- "colortable"
							legend <- FALSE
						}
					} else {
						type <- "colortable"
						legend <- FALSE
					}
				} else if (is.factor(x)) {
					type <- "factor"
					cats <- levels(x)[[1]]
				} else if (is.bool(x)) {
					type <- "factor"
					levels(x) <- data.frame(id=0:1, value=c("False", "True"))
					cats <- cats(x)[[1]]
				} else {
					type <- "depends"
				}
			} else {
				type <- match.arg(type, c("continuous", "classes", "interval", "rgb"))
			}
		}
		if ((type == "classes") && is.null(levels) && is.factor(x)) {
			type <- "factor"
			cats <- cats(x)[[1]]
			if (has.colors(x)) {
				coltab <- coltab(x)[[1]]
			}
		}

		x <- .prep.plot.data(x, type=type, cols=col, mar=mar, draw=TRUE, plg=plg, pax=pax, legend=isTRUE(legend), axes=isTRUE(axes), coltab=coltab, cats=cats, interpolate=smooth, levels=levels, range=range, colNA=colNA, alpha=alpha, reset=reset, grid=grid, sort=sort, decreasing=decreasing, ext=ext, all_levels=all_levels, breaks=breaks, breakby=breakby, add=add, buffer=buffer, background=background, box=box, maxcell=maxcell, clip=clip, ...)

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
	function(x, y, main, mar=NULL, nc, nr, maxnl=16, maxcell=500000, add=FALSE, ...)  {

		if (has.RGB(x)) {
			if (missing(main)) main = ""
			p <- plot(x, -1, main=main, mar=mar, maxcell=maxcell, add=add, ...)
			return(invisible(p))
		}

		nl <- max(1, min(nlyr(x), maxnl))

		if (add && (nl > 1)) {
			nl <- 1
			warn("plot", "adding the first layer of x")
		}

		if (nl==1) {
			if (missing(main)) main = ""
			out <- plot(x, 1, maxcell=maxcell, main=main[1], mar=mar, add=add, ...)
			return(invisible(out))
		}

		nrnc <- .get_nrnc(nr, nc, nl)
		old.par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old.par))
		if (is.null(mar)) {
			mar=c(1.5, 1, 2.5, 3)
		}
		graphics::par(mfrow=nrnc)
		maxcell= 2 * maxcell / nl

		if (missing("main")) {
			tm <- time(x)
			if ((!any(is.na(tm))) && (length(unique(tm)) > 1)) {
				main <- as.character(time(x))
			} else {
				main <- names(x)
			}
		} else {
			main <- rep_len(main, nl)
		}
		for (i in 1:nl) {
			plot(x, i, main=main[i], mar=mar, maxcell=maxcell, add=add, ...)
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


setMethod("plotRGB", signature(x="SpatRaster"),
	function(x, r=1, g=2, b=3, a=NULL, scale=NULL, mar=0, stretch=NULL, smooth=TRUE, 
		colNA="white", alpha=NULL, bgalpha=NULL, zlim=NULL, zcol=FALSE, axes=FALSE, ...) {
	
	x <- x[[c(r, g, b, a)]]
	RGB(x) <- 1:nlyr(x)
	plot(x, -1, scale=scale, mar=mar, stretch=stretch, smooth=smooth, colNA=colNA,
		alpha=alpha, bgalpha=bgalpha, zlim=zlim, zcol=zcol, axes=axes, ...)
}
)

