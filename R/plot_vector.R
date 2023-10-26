# Author: Robert J. Hijmans
# Date :  June 2019
# Version 1.0
# License GPL v3


setMethod("dots", signature(x="SpatVector"),
	function(x, field, size,  ...) {
		reset.clip()

		n <- length(x)
		if (n < 1) return(NULL)
		#method <- match.arg(tolower(method), c("regular", "random"))
		if (is.character(field)) {
			stopifnot(field %in% names(x))
		} else {
			stopifnot(field > 0 && field <= ncol(x))
		}
		stopifnot(is.numeric(x[[field,drop=TRUE]]))
		field <- x[[field,drop=TRUE]]
		size <- size[1]
		stopifnot(size > 0)
		d <- round(field / size)
		d[d < 1 | is.na(d)] <- 0
		i <- d > 0;
		if (sum(i) == 0) {
			error("dots", "'size' is too small")
		}
		s <- spatSample(x[i], d[i], method="random")
		if (.Device  != "null device") {
			try(points(s, ...), silent=TRUE)
		}
		invisible(s)
	}
)



.plotLines <- function(x, out, lty=1, lwd=1, ...) {

	n <- nrow(x)
	if (n == 0) return(out)
#	cols <- out$cols
#	if (is.null(cols)) cols = rep("black", n)

#	g <- lapply(x@cpp$linesList(), function(i) { names(i)=c("x", "y"); i } )

#	g <- geom(x, df=TRUE)
#	g <- split(g, g[,1])
#	g <- lapply(g, function(x) split(x[,3:4], x[,2]))
#	n <- length(g)

	g <- x@cpp$linesList()
	lty <- rep_len(lty, n)
	lwd <- rep_len(lwd, n)
	for (i in 1:n) {
		if (!is.null(g[[i]])) {
			names(g[[i]]) = c("x", "y")
			graphics::plot.xy(g[[i]], type="l", lty=lty[i], col=out$main_cols[i], lwd=lwd[i], ...)
		}
	}

#	for (i in 1:n) {
#		for (j in 1:length(g[[i]])) {
#			lines(g[[i]][[j]], col=out$main_cols[i], lwd=lwd[i], lty=lty[i], ...)
#		}
#	}


	out$leg$lwd <- lwd
	out$leg$lty <- lty
	out
}

.plotPolygons <- function(x, out, lty=1, lwd=1, density=NULL, angle=45, ...) {

	n <- nrow(x)
	if (n == 0) return(out)
	if ((length(out$main_cols) == 0) && (length(out$leg$border) <= 1) && (length(lty) <= 1) && (length(lwd) <= 1) && (is.null(density))) {
		if (is.null(out$leg$border)) out$leg$border <- "black"
		lines(x, col=out$leg$border, lty=lty, lwd=lwd, ...)
		return(out)
	}

#		cols <- .getCols(length(x), col, alpha)
#		out <- list(main_cols=cols)

	if (!is.null(out$leg$border)) {
		out$leg$border <- rep_len(out$leg$border, n)
	} else {
		out$leg$border <- NA
	}
	if (!is.null(density)) {
		out$leg$density <- rep_len(density, n)
		out$leg$angle <- rep_len(angle, n)
	}
	if (!is.null(lty)) out$leg$lty <- rep_len(lty, n)
	if (!is.null(lwd)) out$leg$lwd <- rep_len(lwd, n)

	if (!is.null(out$main_cols)) {
		out$main_cols <- rep_len(out$main_cols, n)
	}

#	g <- geom(x, df=TRUE)
#	g <- split(g, g[,1])
#	g <- lapply(g, function(y) split(y[,3:5], y[,2]))
#	w <- getOption("warn")
#	on.exit(options("warn" = w))
#	for (i in 1:length(g)) {
#		gg <- g[[i]]
#		for (j in 1:length(gg)) {
#			a <- gg[[j]]
#			if (any(is.na(a))) next
#			if (any(a[,3] > 0)) {
#				a <- split(a[,1:2,drop=FALSE], a[,3])
#				a <- lapply(a, function(i) rbind(i, NA))
#				a <- do.call(rbind, a )
#				a <- a[-nrow(a), ]
#				# g[[i]][[1]] <- a
#			}

	g <- x@cpp$polygonsList()
	if (is.null(out$leg$density)) {
		for (i in seq_along(g)) {
			for (j in seq_along(g[[i]])) {
				if (any(is.na(g[[i]][[j]]))) next
				graphics::polypath(g[[i]][[j]][[1]], g[[i]][[j]][[2]], col=out$main_cols[i], rule = "evenodd", border=out$leg$border[i], lwd=out$leg$lwd[i], lty=out$leg$lty[i], ...)
			}
		}
	} else {
		for (i in 1:length(g)) {
			for (j in seq_along(g[[i]])) {
				if (any(is.na(g[[i]][[j]]))) next
				graphics::polygon(g[[i]][[j]][[1]], g[[i]][[j]][[2]], col=out$main_cols[i], density=out$leg$density[i], angle=out$leg$angle[i], border=NA, lwd=out$leg$lwd[i], lty=out$leg$lty[i], ...)
				graphics::polypath(g[[i]][[j]][[1]], g[[i]][[j]][[2]], col=NA, rule="evenodd", border=out$leg$border[i], lwd=out$leg$lwd[i], lty=out$leg$lty[i], ...)
			}
		}
	}
	invisible(out)
}


.vplot <- function(x, out, xlab="", ylab="", pch=16, lty=1, lwd=1, ...) {

	if (out$leg$geomtype == "points") {
		points(x, col=out$main_cols, cex=out$cex, pch=pch, ...)
		#if (!out$add) {
		#	e <- out$lim
		#}
		out$leg$pch = pch
		out$leg$pt.cex = out$cex
	} else if (out$leg$geomtype == "polygons") {
		out <- .plotPolygons(x, out, density=out$leg$density, angle=out$leg$angle, lty=lty, lwd=lwd, ...)
	} else {
		# out <- .plotLines(x, out, ...)
		lines(x, col=out$main_cols, lty=lty, lwd=lwd, ...)
		out$leg$lwd = lwd
		out$leg$lty = lty
	}
	out
}


.getCols <- function(n, cols, alpha=NULL) {
	if (is.null(cols)) { 
		return(cols)
	}
	if (inherits(cols, "function")) {
		cols <- cols(n)
	} else {
		ncols <- length(cols)
		if (ncols > n) {
			steps <- ncols/n
			i <- round(seq(1, ncols, steps))
			cols <- cols[i]
		} else if (ncols < n) {
			cols <- rep_len(cols, n)
		}
	}
	if (!is.null(alpha)) {
		if (alpha[1] < 1 && alpha[1] >= 0) {
			cols <- grDevices::rgb(t(grDevices::col2rgb(cols)), alpha=alpha[1]*255, maxColorValue=255)
		}
	}
	cols
}

.vect.legend.none <- function(out) {
	#if (out$leg$geomtype == "points") {
	#	out$main_cols <- .getCols(out$ngeom, out$cols, 1)
	#} else {
	#	out$cols <- .getCols(out$ngeom, out$cols)
	#}
	out$main_cols <- out$cols
	out
}

.vect.legend.classes <- function(out) {

	if (isTRUE(out$legend_sort)) {
		out$uv <- sort(out$uv, decreasing=out$legend_sort_decreasing)
	} else {
		out$uv <- out$uv[!is.na(out$uv)]
	}
	ucols <- .getCols(length(out$uv), out$cols, out$alpha)

	i <- match(out$v, out$uv)
	out$cols <- ucols
	out$main_cols <- ucols[i]

	if (!is.null(out$colNA)) {
		out$main_cols[is.na(out$main_cols)] <- out$colNA
	}

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

	if (is.null(out$leg$x)) { # && is.null(out$leg$ext)) {
		out$leg$x <- "default"
	}

	out
}


.vect.legend.continuous <- function(out) {

	z <- stats::na.omit(out$v)
	n <- length(z)
	if (n == 0) error("plot", "no values")
	if (!is.numeric(out$v)) {
		out$v <- as.integer(as.factor(out$v))
		z <- stats::na.omit(out$v)
		n <- length(z)
	}
	#out$range <- range(z)

	interval <- (out$range[2]-out$range[1])/(length(out$cols)-1)
	breaks <- out$range[1] + interval * (0:(length(out$cols)-1))

	out$legend_type <- "continuous"
	if (is.null(out$levels)) {
		out$levels <- 5
	}
	if (is.null(out$leg$digits)) {
		dif <- diff(out$range)
		if (dif == 0) {
			out$leg$digits <- 0;
		} else {
			out$leg$digits <- max(0, -floor(log10(dif/10)))
		}
	}

	if (is.null(out$leg$loc)) out$leg$loc <- "right"

	brks <- seq(out$range[1], out$range[2], length.out = length(out$cols))
	grps <- cut(out$v, breaks = brks, include.lowest = TRUE)
	out$main_cols <- out$cols[grps]

	out
}


.vect.legend.interval <- function(out, dig.lab=3, ...) {

	nmx <- length(out$uv)
	if (!is.numeric(out$v)) {
		out$v <- as.integer(as.factor(out$v))
	}

	if (is.null(out$breaks)) {
		out$breaks <- min(5, nmx)
	}

	if (length(out$breaks) == 1) {
		out$breaks <- .get_breaks(out$v, out$breaks, out$breakby, out$range)
	}

	fz <- cut(out$v, out$breaks, include.lowest=TRUE, right=FALSE, dig.lab=dig.lab)
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
		if (length(out$leg$legend) != nlevs) {
			warn("plot", "legend does not match number of levels")
			out$leg$legend <- rep_len(out$leg$legend, nlevs)
		}
	} else {
		levs <- gsub("]", "", gsub(")", "", gsub("\\[", "", levs)))
		levs <- paste(levs, collapse=",")
		m <- matrix(as.numeric(unlist(strsplit(levs, ","))), ncol=2, byrow=TRUE)
		m <- apply(m, 1, function(i) paste(i, collapse=" - "))
		out$leg$legend <- m
	}

	if (is.null(out$leg$x)) { # && is.null(out$leg$ext)) {
		out$leg$x <- "default"
	}
	out$main_cols <- out$cols[out$vcut]
	if (!is.null(out$colNA)) {
		out$main_cols[is.na(out$main_cols)] <- out$colNA
	}
	out
}


.plot.vect.map <- function(x, out, ...) {

	if ((!out$add) & (!out$legend_only)) {
		if (!any(is.na(out$mar))) { graphics::par(mar=out$mar) }
		plot(out$lim[1:2], out$lim[3:4], type="n", xlab="", ylab="", asp=out$asp, xaxs="i", yaxs="i", axes=FALSE, main="")
		if (!is.null(out$background)) {
			graphics::rect(out$lim[1], out$lim[3], out$lim[2], out$lim[4], col=out$background)
		}
	}
	if (isTRUE(out$blank)) return(out)
	
	nuq <- length(out$uv)
	if (out$legend_type == "none") {
		out <- .vect.legend.none(out)
	} else if (out$legend_type == "classes") {
		out <- .vect.legend.classes(out)
	} else if (out$legend_type == "interval") {
		if (nuq < 2) {
			out <- .vect.legend.classes(out, ...)
		} else {
			out <- .vect.legend.interval(out, dig.lab=out$dig.lab)
		}
	} else if (out$legend_type == "depends") {
		if (!is.null(out$breaks)) {
			out <- .vect.legend.interval(out, dig.lab=out$dig.lab)			
		} else if (nuq < 11) {
			out <- .vect.legend.classes(out)
		} else if (!is.numeric(out$uv)) {
			#if (nuq < 21)
			out <- .vect.legend.classes(out)
		} else {
			out <- .vect.legend.interval(out, dig.lab=out$dig.lab)
		}
	} else {
		if (nuq == 1) {
			out <- .vect.legend.classes(out)
		} else {
			out <- .vect.legend.continuous(out)
			out$leg$density <- NULL
		}
	}
	if (!out$legend_only) {
		if (!out$add) {
			try(set.clip(out$lim, out$lonlat))
		}
		out <- .vplot(x, out, ...)
	}

	if (out$axes) {
		out <- .plot.axes(out)
	}

	if (out$legend_draw) {
		if (out$legend_type == "continuous") {
			out$legpars <- do.call(.plot.cont.legend, list(x=out))
		} else {
			if (out$add) {
				if (out$clip) {
					out$leg$plotlim <- unlist(get.clip()[1:4])
				} else {
					out$leg$plotlim <- graphics::par("usr")
				}
				if (is.null(out$leg$plotlim)) {
					out$leg$plotlim <- out$lim			
				}
			} else {
				if (out$clip) {
					out$leg$plotlim <- out$lim
				} else {
					out$leg$plotlim <- graphics::par("usr")
				}
			}			
			out$legpars <- do.call(.plot.class.legend, out$leg)
		}
	}
	if (isTRUE(out$box)) { 
		if (out$clip) {
			lines(ext(out$lim))	
		} else {
			lines(ext(graphics::par("usr")))		
		}
	}

	if (out$main != "") {
		posx <- out$lim[1] + diff(out$lim[1:2])/2
		text(posx, out$lim[4], out$main, pos=3, offset=out$line.main, cex=out$cex.main, 
			font=out$font.main, col=out$col.main, xpd=TRUE)
	}

	if (!out$add) {
		try(set.clip(out$lim, out$lonlat))
	}
	out
}


.prep.vect.data <- function(x, y, type=NULL, cols=NULL, mar=NULL, legend=TRUE,
	legend.only=FALSE, levels=NULL, add=FALSE, range=NULL, breaks=NULL, breakby="eqint",
	xlim=NULL, ylim=NULL, colNA=NA, alpha=NULL, axes=TRUE, buffer=TRUE, background=NULL,
	pax=list(), plg=list(), ext=NULL, grid=FALSE, las=0, sort=TRUE, decreasing=FALSE, values=NULL,
	box=TRUE, xlab="", ylab="", cex.lab=0.8, line.lab=1.5, yaxs="i", xaxs="i", main="", cex.main=1.2, line.main=0.5, font.main=graphics::par()$font.main, col.main = graphics::par()$col.main, 
	density=NULL, angle=45, border="black", dig.lab=3, cex=1, clip=TRUE, leg_i=1, asp=NULL, ...) {

	out <- list()
	out$blank <- FALSE
	if ((y == "") && (is.null(values))) {
		if (is.null(type)) type <- "none"
		if (type == "n") {
			out$blank <- TRUE
		}
		type <- "none"
		plg <- list()
	}
	
	e <- as.vector(ext(x))
	if (any(is.na(e))) {
		error("plot", "SpatVector has no valid geometries")
	}
	if (e[1] == e[2]) {
		e[1] = e[1] - 0.5
		e[2] = e[2] + 0.5
	}
	if (e[3] == e[4]) {
		e[3] = e[3] - 0.5
		e[4] = e[4] + 0.5
	}
	
	out$lim <- out$ext <- e
	if ((!is.null(ext)) || (!is.null(xlim)) || (!is.null(ylim))) {
		if (!is.null(ext)) {
			ext <- ext(ext)
			if (!is.null(xlim) | !is.null(ylim)) {
				x <- crop(x, ext)
			} else {
				x <- x[ext, ]
			}
			out$ext <- as.vector(ext(x))
			out$lim <- as.vector(ext)
		}
		if (!is.null(xlim)) {
			stopifnot(length(xlim) == 2)
			out$lim[1:2] <- sort(xlim)
		}
		if (!is.null(ylim)) {
			stopifnot(length(ylim) == 2)
			out$lim[3:4] <- sort(ylim)
		}
	} 
	out$ngeom <- nrow(x)
	out$clip <- isTRUE(clip)

	if (buffer) {
		dx <- diff(out$lim[1:2]) / 50
		dy <- diff(out$lim[3:4]) / 50
		out$lim[1:2] <- out$lim[1:2] + c(-dx, dx)
		out$lim[3:4] <- out$lim[3:4] + c(-dy, dy)
	}
	
	out$main <- main
	if ((!is.expression(main)) && (is.null(out$main) || any(is.na(out$main)))) out$main <- ""
	out$cex.main  <- cex.main
	out$font.main <- font.main
	out$col.main  <- col.main
	out$line.main <- line.main
	out$dig.lab <- dig.lab

	out$box <- isTRUE(box)
	out$add <- isTRUE(add)
	out$axes <- isTRUE(axes)
	out$xlab <- xlab
	out$ylab <- ylab
	out$axs <- as.list(pax)
	out$cex <- cex
	if (is.null(out$axs$las)) out$axs$las <- las
	if (is.null(out$axs$cex.lab)) out$axs$cex.lab <- cex.lab
	if (is.null(out$axs$line.lab)) out$axs$line.lab <- line.lab
	
	out$draw_grid <- isTRUE(grid)
	out$leg <- as.list(plg)
	out$leg$geomtype <- geomtype(x)
	
	out$leg$density <- density
	out$leg$angle <- angle
	out$leg$border <- border
	
	
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
	
	out$breaks <- breaks
	out$breakby <- breakby
	out$background <- background
	if (is.null(values)) {
		v <- unlist(x[, y, drop=TRUE], use.names=FALSE)
	} else {
		if (inherits(values, "data.frame")) {
			if (ncol(values) == 2) {
				xname = names(values)[1]
				if (xname %in% names(x)) {
					i <- match(x[[xname,drop=TRUE]], values[[1]])
					v <- values[[2]][i]
				} else {
					error("plot", paste(xname, "is not a name in x"))
				}
			} else {
				values <- values[[1]]
			}
		} else {
			v <- as.vector(values)
		}
		v <- rep_len(v, nrow(x))
	}
	if (is.factor(v)) v <- as.character(v)
	if (is.numeric(v)) {
		v[!is.finite(v)] <- NA
	}
	if (!is.null(range)) {
		range <- sort(range)
		v[v < range[1]] <- NA
		v[v > range[2]] <- NA
		if (all(is.na(v))) {
			v <- NULL
			y <- ""
			type = "none"
		} else {
			out$range <- range
		}
		out$range_set <- TRUE
	} else {
		if (!is.null(v)) {
			out$range <- range(v, na.rm=TRUE)
		}
		out$range_set <- FALSE
	}
	out$v <- v

	if (!is.logical(sort)) {
		out$uv <- unique(sort)
		out$legend_sort <- FALSE
	} else {
		out$uv <- unique(out$v)
		out$legend_sort <- isTRUE(sort)
		out$legend_sort_decreasing <- isTRUE(decreasing)
	}

	if (is.null(type)) {
		type <- "depends"
	} else {
		type <- match.arg(type, c("continuous", "classes", "interval", "depends", "none"))
	}
	out$levels <- levels
	
	if (type=="none") {
		legend <- FALSE
		legend_only <- FALSE
	}
	out$legend_type <- type

	if (is.null(cols)) {
		if (type == "none") {
			if (out$leg$geomtype %in% c("lines", "points")) {
				cols <- "black"
			}
		} else {
			cols <- rev(grDevices::rainbow(100, start=.1, end=0.9))
		}
	}
	if (!is.null(alpha)) {
		alpha <- clamp(alpha[1], 0, 1)
		cols <- grDevices::rgb(t(grDevices::col2rgb(cols)), alpha=alpha*255, maxColorValue=255)
	} else {
		alpha <- 1
	}

	out$alpha <- alpha
	out$cols <- cols
	out$legend_draw <- isTRUE(legend)
	out$legend_only <- isTRUE(legend.only)
	out$leg$leg_i <- leg_i


	if (is.null(mar)) {
		if (out$legend_draw) {
			mar=c(3.1, 3.1, 2.1, 7.1)
		} else {
			mar=c(3.1, 3.1, 2.1, 2.1)
		}
	}
	out$mar <- rep_len(mar, 4)

	out$skipNA <- TRUE
	if (!is.null(colNA)) {
		if (!is.na(colNA)) {
			out$colNA <- grDevices::rgb(t(grDevices::col2rgb(colNA)), alpha=alpha*255, maxColorValue=255)
			out$r[is.na(out$r)] <- out$colNA
			out$skipNA <- FALSE
 		} else {
			out$colNA <- NULL
		}
	}

	.plot.vect.map(x, out, ...)
}


setMethod("plot", signature(x="SpatVector", y="character"),
	function(x, y, col=NULL, type=NULL, mar=NULL, add=FALSE, legend=TRUE, axes=!add,
	main, buffer=TRUE, background=NULL, grid=FALSE, ext=NULL, 
	sort=TRUE, decreasing=FALSE, plg=list(), pax=list(), nr, nc, colNA=NA, 
	alpha=NULL, box=axes, clip=TRUE, ...) {

		old.mar <- graphics::par()$mar
		on.exit(graphics::par(mar=old.mar))

		if (nrow(x) == 0) {
			error("plot", "SpatVector has zero geometries")
		}
		if (add) reset.clip()

		y <- trimws(y)
		if (any(is.na(match(y, c("", names(x)))))) {
			i <- is.na(match(y, names(x)))
			error("plot", paste(paste(y[i], collapse=",")), " is not a name in x")
		}
		nrnc <- c(1,1)
		if (length(y) > 1) {
			nrnc <- .get_nrnc(nr, nc, length(y))
			old.par <- graphics::par(no.readonly =TRUE)
			on.exit(graphics::par(old.par))
			graphics::par(mfrow=nrnc)
		}
		if (is.character(legend)) {
			plg$x <- legend
			legend <- TRUE
		}


		for (i in 1:length(y)) {
			if (length(y) > 1) {
				if (missing("main")) {
					main <- y
				} else {
					main <- rep_len(main, length(y))
				}
				newrow <- (nrnc[2] == 1) | ((i %% nrnc[2]) == 1)
				lastrow <- i > (prod(nrnc) - nrnc[2])
				if (lastrow) {
					if (newrow) {
						pax$side <- 1:2
					} else {
						pax$side <- 1
					}
				} else if (newrow) {
					pax$side <- 2
				} else {
					pax$side <- 0
				}
			} else if (missing("main")) {
				main <- ""
			}

			if (missing(col)) col <- NULL

			out <- .prep.vect.data(x, y[i], type=type, cols=col, mar=mar, plg=plg, pax=pax, legend=isTRUE(legend), add=add, axes=axes, main=main[i], buffer=buffer, background=background, grid=grid, ext=ext, sort=sort, decreasing=decreasing, colNA=colNA, alpha=alpha, box=box, clip=clip, leg_i=i, ...)
		}
		invisible(out)
	}
)


setMethod("plot", signature(x="SpatVector", y="numeric"),
	function(x, y, ...)  {
		y <- round(y)
		if (any(y > ncol(x))) {
			error("plot", paste("x only has", ncol(x), " columns"))
		}
		y[y<0] <- 0
		y <- c("", names(x))[y+1]
		out <- plot(x, y, ...)
		invisible(out)
	}
)


#setMethod("plot", signature(x="SpatVector", y="data.frame"),
#	function(x, y, values=NULL, ...)  {
#		out <- plot(x, "", values=y, ...)
#		invisible(out)
#	}
#)

setMethod("plot", signature(x="SpatVector", y="missing"),
	function(x, y, values=NULL, ...)  {
		invisible( plot(x, "", values=values, ...) )
	}
)


setMethod("plot", signature(x="SpatVectorProxy", y="missing"),
	function(x, y, ...)  {
		plot(ext(x), ...)
	}
)

