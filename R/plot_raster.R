

.as.raster.continuous <- function(out, x, type) {

	Z <- as.matrix(x, wide=TRUE)
	Z[is.nan(Z) | is.infinite(Z)] <- NA
	Z[] <- round(Z, 12)
	
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

	uzi <- unique(z)

	if (type == "depends") {
		if (length(uzi) < 9) {
			return (.as.raster.classes(out, x, Z=Z))
		}
	} else if ((length(uzi) == 1) && is.null(out$range)) {
		return (.as.raster.classes(out, x, Z=Z))
	}

	if (is.null(out$range)) {
		out$range <- range(z)
#		out$fill_range <- FALSE
	} else {
		stopifnot(length(out$range) == 2)
		if (out$fill_range) {
			out$range_filled <- c(FALSE, FALSE)
			if (!is.na(out$range[1])) {
				if (out$range[1] > min(z)) {
					out$range_filled[1] <- TRUE
					Z[ Z < out$range[1] ] <- out$range[1]
				} 
			} else {
				out$range[1] <- min(z, na.rm=TRUE)
			}
			if (!is.na(out$range[2])) {
				if (out$range[2] < max(z)) {
					Z[ Z > out$range[2] ] <- out$range[2]
					out$range_filled[2] <- TRUE
				} 
			} else {
				out$range[2] <- max(z, na.rm=TRUE)
			}
		} else {
			if (all(is.na(out$range))) {
				out$range <- range(z)
			} else if (is.na(out$range[1])) {
				out$range[1] <- min(z)
			} else if (is.na(out$range[2])) {
				out$range[2] <- max(z)
			}
		}
		
		if (!any(out$range_filled)) out$fill_range <- FALSE
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

	if (length(breaks) == 1) {
		Z[] <- out$cols[ceiling(length(out$cols)/2)]
	} else {
		Z[] <- out$cols[as.integer(cut(as.numeric(Z), breaks, include.lowest=TRUE, right=FALSE))]
	}
	
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


.as.raster.classes <- function(out, x, Z=NULL, ...) {

	if (is.null(Z)) {
		Z <- as.matrix(x, wide=TRUE)
		Z[is.nan(Z) | is.infinite(Z)] <- NA
	}
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
		#i <- match(Z, as.numeric(levs))
		i <- match(Z, as.numeric(out$cols[,1]))
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
	if (!is.null(out$leg$order)) {
		ord <- match(out$leg$order, out$leg$legend)
		out$leg$legend <- out$leg$legend[ord]
		out$leg$fill <- out$leg$fill[ord]
	} else if (isTRUE(out$leg$sort)) {
		ord <- order(out$leg$legend) #, decreasing=out$leg$reverse)
		out$leg$legend <- out$leg$legend[ord]
		out$leg$fill <- out$leg$fill[ord]
	}
	out
}



.as.raster.interval <- function(out, x, ...) {

	out$legend_type <- "classes"

	if (NCOL(out$cols) == 3) {
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
		nZ <- length(unique(na.omit(as.vector(Z))))	
		if (nZ <= 1) {
			return(.as.raster.classes(out, x, Z=Z, ...))
		}

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
		x$zebra <- FALSE
	} else if (!x$legend_only) {
		old.mar <- graphics::par()$mar
		if (!any(is.na(x$mar))) { graphics::par(mar=x$mar) }
		if (x$reset) on.exit(graphics::par(mar=old.mar))
		
		if (x$zebra) {
			width <- rep(min(diff(x$lim[1:2]), diff(x$lim[3:4])) / 100, 2) * x$zebra.cex
			if (x$lonlat) {
				asp <- 1/cos((mean(x$lim[3:4]) * pi)/180)
				width[2] <- width[2] / asp
			}
			x$lim[1:2] <- x$lim[1:2] + c(-width[1], width[1])
			x$lim[3:4] <- x$lim[3:4] + c(-width[2], width[2])
		}		
		
		arglist <- c(list(x=x$lim[1:2], y=x$lim[3:4], type="n", xlab="", ylab="", asp=x$asp, xaxs=x$xaxs, yaxs=x$yaxs, axes=FALSE), x$dots)
		do.call(plot, arglist)
		if (!is.null(x$background)) {
			graphics::rect(x$lim[1], x$lim[3], x$lim[2], x$lim[4], col=x$background, border=x$box)			
		}
	}
	try(set.clip(x$lim, x$lonlat))
	if (!x$values) {
		#if (!x$add) try(set.clip(x$lim, x$lonlat))
		return(x)
	}
	if (!x$legend_only) {
		x <- .plot.axes(x)
		graphics::rasterImage(x$r, x$ext[1], x$ext[3], x$ext[2], x$ext[4], angle = 0, interpolate = x$interpolate)
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

	if (x$zebra) {
		try(set.clip(x$lim, x$lonlat))
		zebra(width=width, x=x$axs$xat, y=x$axs$yat, col=x$zebra.col)
		x$lim[1:2] <- x$lim[1:2] + c(width[1], -width[1])
		x$lim[3:4] <- x$lim[3:4] + c(width[2], -width[2])
	}		

	if (isTRUE(x$box)) { 
		if (x$clip) {
			lines(ext(x$lim))	
		} else {
			lines(ext(graphics::par("usr")))		
		}
	}
	
	if ((!x$add) && (!x$legend_only)) {
		plot_main(x)
	}
	if (!x$add) {
		try(set.clip(x$lim, x$lonlat))
	} else {
		reset.clip()
	}	
	invisible(x)
}



.prep.plot.data <- function(x, type, cols, mar=NULL, draw=FALSE,
  interpolate=FALSE, legend=TRUE, legend.only=FALSE, pax=list(), plg=list(),
  levels=NULL, add=FALSE, range=NULL, fill_range=FALSE, breaks=NULL, breakby="eqint",
  coltab=NULL, cats=NULL, xlim=NULL, ylim=NULL, ext=NULL, colNA=NA, alpha=NULL, reset=FALSE,
  sort=TRUE, reverse=FALSE, grid=FALSE, las=0, all_levels=FALSE, decimals=NULL, background=NULL,
  xlab="", ylab="", cex.lab=0.8, line.lab=1.5, asp=NULL, yaxs="i", xaxs="i", 
  main="", cex.main=1.2, line.main=0.5, font.main=graphics::par()$font.main, col.main = graphics::par()$col.main, loc.main=NULL, 
  sub = "", font.sub=1, cex.sub=0.8*cex.main, line.sub =1.75,  col.sub=col.main, loc.sub=NULL,
  halo=FALSE, hc="white", hw=0.1, axes=TRUE, box=TRUE, zebra=FALSE, zebra.cex=1, zebra.col=c("black", "white"), 
  maxcell=500000, buffer=FALSE, clip=TRUE, 
  # for rgb 
  stretch=NULL, scale=NULL, bgalpha=NULL, zlim=NULL, zcol=NULL, overview=NULL, 
#catch and kill
  cex=1, decreasing=FALSE, font=NULL,
  ...) {


	# backwards compatibility
	reverse <- reverse | decreasing
	out <- list()
	e <- out$lim <- out$ext <- as.vector(ext(x))
	hadWin <- hasWin <- FALSE
	if (add && is.null(ext)) {
		ext <- unlist(get.clip())[1:4]
	}
	
	if ((!is.null(ext)) || (!is.null(xlim)) || (!is.null(ylim))) {
		if (!is.null(ext)) {
			ext <- ext(ext)
			#e <- as.vector(align(ext, x))
			e <- as.vector(ext)
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
			
		if (is.null(overview)) {	
			if (grepl("https://", tolower(sources(x))[1])) {
				overview <- TRUE
			} else {
				overview <- FALSE
			}
		}
	
		if (inherits(alpha, "SpatRaster")) {
			if (nlyr(alpha) > 1) {
				alpha <- alpha[[1]]
			}
#			alpha <- spatSample(alpha, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
			alpha <- sampleRaster(alpha, maxcell, method="regular", replace=FALSE, ext=NULL, warn=FALSE, overview=overview)
		}
#		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
		x <- sampleRaster(x, maxcell, method="regular", replace=FALSE, ext=NULL, warn=FALSE, overview=overview)

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
	out$cols <- cols
	out$rgb$stretch <- stretch
	if (is.null(scale)) {	
		out$rgb$scale <- 255
	} else {
		out$rgb$scale <- scale
	}
	if (!is.null(alpha)) {
		if (!inherits(alpha, "SpatRaster")) {
			out$alpha <- alpha[1]
			if ((alpha < 0) || (alpha > 1)) {
				warn("plot", "alpha should be between 0 and 1")
				out$alpha <- 255
			} else {
				out$alpha <- out$alpha[1] * 255
			}
			out <- hexcols(out)
		}
	} else {
		out$alpha <- 255
		if (!is.null(scale)) {
			out$alpha <- out$alpha * scale / 255
		}
		out <- hexcols(out)
	}
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
	out$halo.main.hc <- hc
	out$halo.main.hw <- hw

	out$loc.main  <- loc.main
	out$cex.main  <- cex.main
	out$font.main <- font.main
	out$col.main  <- col.main
	out$line.main <- line.main

	out$sub <- sub
	out$loc.sub <- loc.sub
	out$cex.sub <- cex.sub
	out$font.sub <- font.sub
	out$col.sub <- col.sub
	out$line.sub <- line.sub

	out$axes <- axes
	out$xaxs <- xaxs
	out$yaxs <- yaxs
	out$xlab <- xlab
	out$ylab <- ylab 
	out$coltab <- coltab
	out$cats <- cats
	out$breaks <- breaks
	out$breakby <- breakby
	out$interpolate <- FALSE
	out$background <- background
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)
	if (!is.logical(sort)) {
		out$leg$order <- sort
		out$leg$sort <- FALSE
	} else {
		out$leg$sort <- isTRUE(sort)
	}
	out$leg$reverse <- isTRUE(reverse)
	out$box <- isTRUE(box)
	out$zebra <- isTRUE(zebra)
	out$zebra.cex <- zebra.cex
	out$zebra.col <- zebra.col
	
#	if (!is.null(out$leg$loc)) {
#		out$leg$x <- out$leg$loc
#		out$leg$loc <- NULL
#	}


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
			out$fill_range <- fill_range
			out <- .as.raster.continuous(out, x, type)
		}

		out$mar <- mar
		out <- get_mar(out)
		
		if (!is.null(colNA)) {
			if (!is.na(colNA) && out$values) {
				out$colNA <- grDevices::rgb(t(grDevices::col2rgb(colNA)), alpha=out$alpha, maxColorValue=255)
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
	function(x, y=1, col, type=NULL, mar=NULL, legend=TRUE, axes=!add, plg=list(), pax=list(), maxcell=500000, smooth=FALSE, range=NULL, fill_range=FALSE, levels=NULL, all_levels=FALSE, breaks=NULL, breakby="eqint", fun=NULL, colNA=NULL, alpha=NULL, sort=FALSE, reverse=FALSE, grid=FALSE, zebra=FALSE, ext=NULL, reset=FALSE, add=FALSE, buffer=FALSE, background=NULL, box=axes, clip=TRUE, overview=NULL, ...) {

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
			plot(x, col=col, type=type, mar=mar, legend=legend, axes=axes, plg=plg, pax=pax, maxcell=2*maxcell/length(y), smooth=smooth, range=range, fill_range=fill_range, levels=levels, all_levels=all_levels, breaks=breaks, breakby=breakby, fun=fun, colNA=colNA, alpha=alpha, zebra=zebra, grid=grid, sort=sort, reverse=reverse, ext=ext, reset=reset, add=add, buffer=buffer, background=background, box=box, clip=clip, overview=overview, ...)
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
			col <- .default.pal()
			#col <- rev(grDevices::terrain.colors(255))
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

		x <- .prep.plot.data(x, type=type, cols=col, mar=mar, draw=TRUE, plg=plg, pax=pax, legend=isTRUE(legend), axes=isTRUE(axes), coltab=coltab, cats=cats, interpolate=smooth, levels=levels, range=range, fill_range=fill_range, colNA=colNA, alpha=alpha, reset=reset, grid=grid, zebra=zebra, sort=sort, reverse=reverse, ext=ext, all_levels=all_levels, breaks=breaks, breakby=breakby, add=add, buffer=buffer, background=background, box=box, maxcell=maxcell, clip=clip, overview=overview, ...)

		add_more(fun, y)
		
		invisible(x)
	}
)


setMethod("plot", signature(x="SpatRaster", y="missing"),
	function(x, y, main, mar=NULL, nc, nr, maxnl=16, maxcell=500000, add=FALSE, plg=list(), pax=list(), ...)  {

		if (has.RGB(x)) {
			if (missing(main)) main = ""
			p <- plot(x, -1, main=main, mar=mar, maxcell=maxcell, add=add, plg=plg, pax=pax, ...)
			return(invisible(p))
		}

		nl <- max(1, min(nlyr(x), maxnl))

		if (add && (nl > 1)) {
			nl <- 1
			warn("plot", "adding the first layer of x")
		}

		if (nl==1) {
			if (missing(main)) main = ""
			out <- plot(x, 1, maxcell=maxcell, main=main[1], mar=mar, add=add, plg=plg, pax=pax, ...)
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
		plg$title <- suppressWarnings(rep(plg$title, length.out=nl))

		for (i in 1:nl) {
			plg$leg_i <- i
			plot(x, i, main=main[i], mar=mar, maxcell=maxcell, add=add, plg=plg, pax=pax, ...)
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

