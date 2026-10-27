


.one.density <- function(x, maxcells=100000, plot=TRUE, ...) {
	d <- values(x)
	d <- stats::density(stats::na.omit(d))
	if (plot) {
		plot(d, ...)
		return(invisible(d))
	} else {
		return(d)
	}
}


setMethod("density", signature(x="SpatRaster"),
	function(x, maxcells=100000, plot=TRUE, main, ...) {
		x <- spatSample(x, maxcells, method="regular", as.raster=TRUE, warn=FALSE)
		res <- list()
		nl <- nlyr(x)
		if (missing(main)) {
			main=names(x)
		} else {
			main <- rep(main, length.out=nl)
		}

		if (nl==1) {
			res[[1]] <- .one.density(x, plot=plot, maxcells=maxcells, main=main, ...)
		} else {

			if (nl > 16) {
			warn("density", "only the first 16 layers are plotted")
				nl <- 16
				x <- x[[1:16]]
			}
			nc <- ceiling(sqrt(nl))
			nr <- ceiling(nl / nc)
			mfrow <- graphics::par("mfrow")
			spots <- mfrow[1] * mfrow[2]
			if (spots < nl) {
				old.par <- graphics::par(no.readonly = TRUE)
				on.exit(graphics::par(old.par))
				graphics::par(mfrow=c(nr, nc))
			}
			for (i in 1:nlyr(x)) {
				res[[i]] <- .one.density(x[[i]], maxcells=maxcells, main=main[i], plot=plot, ...)
			}
		}
		if (plot) return(invisible(res))
		else return(res)
	}
)


setMethod("persp", signature(x="SpatRaster"),
	function(x, maxcells=100000, ...)  {
		x <- spatSample(x, size=maxcells, method="regular", as.raster=TRUE, warn=FALSE)
		value <- t(as.matrix(x, wide=TRUE)[nrow(x):1,])
		y <- yFromRow(x, nrow(x):1)
		x <- xFromCol(x, 1:ncol(x))
		graphics::persp(x=x, y=y, z=value, ...)
	}
)


.plot.filled.contour <- function(x, maxcells=100000, ...) {

	x <- spatSample(x[[1]], maxcells, method="regular", as.raster=TRUE, warn=FALSE)
	X <- xFromCol(x, 1:ncol(x))
	Y <- yFromRow(x, nrow(x):1)
	Z <- t( matrix( values(x), ncol=ncol(x), byrow=TRUE)[nrow(x):1,] )

 	if (is.null(list(...)$asp)) {
		asp <- ifelse(is.lonlat(x, perhaps=TRUE, warn=FALSE), 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180), 1)
		graphics::filled.contour(x=X, y=Y, z=Z, asp=asp, ...)
	} else {
		graphics::filled.contour(x=X, y=Y, z=Z,...)
	}
}


setMethod("contour", signature(x="SpatRaster"),
	function(x, maxcells=100000, filled=FALSE, ...)  {
		if (filled) {
			.plot.filled.contour(x, maxcells=maxcells, ...)
		} else {
			x <- spatSample(x[[1]], maxcells, method="regular", as.raster=TRUE, warn=FALSE)
			if (is.null(list(...)$asp)) {
				asp <- ifelse(is.lonlat(x, perhaps=TRUE, warn=FALSE), 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180), 1)
				graphics::contour(x=xFromCol(x,1:ncol(x)), y=yFromRow(x, nrow(x):1), z=t(as.matrix(x, wide=TRUE)[nrow(x):1,]), asp=asp, ...)
			} else {
				graphics::contour(x=xFromCol(x,1:ncol(x)), y=yFromRow(x, nrow(x):1), z=t(as.matrix(x, wide=TRUE)[nrow(x):1,]), ...)
			}
		}
	}
)


setMethod("as.contour", signature(x="SpatRaster"),
	function(x, maxcells=100000, ...) {
		x <- spatSample(x[[1]], size=maxcells, method="regular", as.raster=TRUE, warn=FALSE)
		z <- grDevices::contourLines(x=xFromCol(x,1:ncol(x)), y=yFromRow(x, nrow(x):1), z=t(as.matrix(x, wide=TRUE)[nrow(x):1,]), ...)
		y <- lapply(1:length(z), function(i) cbind(z[[i]]$level, i, z[[i]]$x, z[[i]]$y))
		y <- do.call(rbind, y)
		y[] <- as.numeric(y)
		u <- unique(y[,1])
		y[,1] <- match(y[,1], u)
		colnames(y)[3:4] <- c("x", "y")
		vect(y, "lines", atts=data.frame(level=u), crs=crs(x))
	}
)



setMethod("pairs", signature(x="SpatRaster"),
	function(x, hist=TRUE, cor=TRUE, use="pairwise.complete.obs",  maxcells=100000, ...) {

		if (nlyr(x) < 2) {
			error("x must have at least two layers")
		}
		if (nlyr(x) < 2) {
			error("x must have at least two layers")
		}

		panelhist <- function(x,...)	{
			usr <- graphics::par("usr")
			on.exit(graphics::par(usr=usr))
			graphics::par(usr = c(usr[1:2], 0, 1.5) )
			h <- hist(x, plot = FALSE)
			breaks <- h$breaks
			nB <- length(breaks)
			y <- h$counts
			y <- y/max(y)
			graphics::rect(breaks[-nB], 0, breaks[-1], y, col="green")
		}

		panelcor <- function(x, y,...) {
			usr <- graphics::par("usr")
			on.exit(graphics::par(usr=usr))
			graphics::par(usr = c(0, 1, 0, 1))
			r <- abs(stats::cor(x, y, use=use))
			txt <- format(c(r, 0.123456789), digits=2)[1]
			text(0.5, 0.5, txt, cex = max(0.5, r * 2))
		}

		if (hist) {dp <- panelhist} else {dp <- NULL}
		if (cor) {up <- panelcor} else {up <- NULL}

		d <- spatSample(x, maxcells, method="regular", as.raster=FALSE, warn=FALSE)

		dots <- list(...)
		cex <- dots$cex
		main <- dots$main
		if (is.null(cex)) cex <- 0.5
		if (is.null(main)) main <- ""

		graphics::pairs(d, main=main, cex=cex, upper.panel=up, diag.panel=dp)
	}
)

.textbox <- function(x, y=NULL, labels, col="black", hc="white", hw=0.1, ... ) {
	xy <- grDevices::xy.coords(x, y)
	hw <- hw[1]
	n <- nchar(labels)
	x0 <- hw * graphics::strwidth("A")
	y0 <- hw * graphics::strheight("A")
	x1 <- n * x0
	v <- vect(cbind(xy$x + c(0, x1), xy$y), crs="local")
	b <- buffer(v, y0)
}


.halo <- function(x, y=NULL, labels, col="black", hc="white", hw=0.1, ... ) {
	if (!is.null(y)) { # recycle
		xy <- grDevices::xy.coords(cbind(x, y))
	} else {
		xy <- grDevices::xy.coords(x, y)
	}
	hw <- hw[1]
	xo <- hw * graphics::strwidth("A")
	yo <- hw * graphics::strheight("A")
	n <- nchar(labels)
	fact <- 100 * max(1, hw*2)
	
	theta <- seq(pi/4, 2*pi, length.out=hw*fact)
	for (i in theta) {
		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=hc, ... )
	}
	graphics::text(xy$x, xy$y, labels, col=col, ... )
}

halo <- function(x, y=NULL, labels, col="black", hc="white", hw=0.1, ... ) {
	.halo(x=x, y=y, labels=labels, col=col, hc=hc, hw=hw, ...) 
}


setMethod("text", signature(x="SpatRaster"),
	function(x, labels, digits=0, halo=FALSE, hc="white", hw=0.1, jitter=0, xpd=TRUE, ...) {
		if (missing(labels)) {
			labels <- 1
		}
		if (length(labels) != ncell(x)) {
			labels <- labels[1]
			if (is.character(labels)) {
				i <- which(labels == names(x))
				if (i == 0) {
					i <- 1
				}
			}
			x <- x[[labels]]
			labels <- as.data.frame(x)[,1]
			p <- as.points(x, values=TRUE)
		} else {
			p <- as.points(x, values=FALSE, na.rm=FALSE)
		}
		xy <- geom(p)[, c("x", "y"), drop=FALSE]
        if (isTRUE(jitter > 0)) {
			xy <- jitter(xy, factor = jitter)
		}
		if (is.factor(labels)) {
			labels <- substr(as.character(labels), 1, max(1, digits))
		} else if (is.numeric(labels)) {
			labels <- as.character(round(labels, digits=digits) )
		}
		if (length(labels) < nrow(xy)) {
			labels <- rep(labels, nrow(xy))
		}
#		if (!overlap) {
#			xy <- getLabelXY(xy, labels, cex)
#		}
		if (halo && (isTRUE(hw > 0))) {
			.halo(xy[,1], xy[,2], labels, hc=hc, hw=hw, xpd=xpd, ...)
		} else {
			text(xy[,1], xy[,2], labels, xpd=xpd, ...)
		}
		invisible(xy)
	}
)


.central_loc <- function(x, labels="", dots=list(), halo=FALSE, hw=0.1, inside=FALSE) {

	.label_boxes <- function(xy, labels, dots, halo=FALSE, hw=0.1) {
		n <- nrow(xy)
		labels <- rep_len(as.character(labels), n)
		cex <- dots$cex
		if (is.null(cex)) cex <- graphics::par("cex")
		cex <- rep_len(cex, n)
		font <- dots$font
		if (is.null(font)) font <- graphics::par("font")
		family <- dots$family
		if (is.null(family)) family <- graphics::par("family")
		vfont <- dots$vfont
		srt <- dots$srt
		if (is.null(srt)) srt <- 0
		srt <- rep_len(srt, n)

		sw <- function(s, c) {
			graphics::strwidth(s, units="user", cex=c, font=font, family=family, vfont=vfont)
		}
		sh <- function(s, c) {
			graphics::strheight(s, units="user", cex=c, font=font, family=family, vfont=vfont)
		}
		w <- h <- numeric(n)
		for (c in unique(cex)) {
			j <- which(cex == c)
			w[j] <- sw(labels[j], c)
			h[j] <- sh(labels[j], c)
		}
		if (halo && isTRUE(hw > 0)) {
			w <- w + 2 * hw * sw(rep("A", n), cex)
			h <- h + 2 * hw * sh(rep("A", n), cex)
		}

		xx <- xy[, 1]
		yy <- xy[, 2]
		ax <- rep(0.5, n)
		ay <- rep(0.5, n)
		pos <- dots$pos
		if (!is.null(pos)) {
			pos <- rep_len(as.integer(pos), n)
			offset <- dots$offset
			if (is.null(offset)) offset <- 0.5
			ox <- offset * sw(rep("m", n), cex)
			oy <- offset * sh(rep("X", n), cex)
			i <- pos == 1
			ax[i] <- 0.5; ay[i] <- 1; yy[i] <- yy[i] - oy[i]
			i <- pos == 2
			ax[i] <- 1; ay[i] <- 0.5; xx[i] <- xx[i] - ox[i]
			i <- pos == 3
			ax[i] <- 0.5; ay[i] <- 0; yy[i] <- yy[i] + oy[i]
			i <- pos == 4
			ax[i] <- 0; ay[i] <- 0.5; xx[i] <- xx[i] + ox[i]
		} else if (!is.null(dots$adj)) {
			adj <- rep_len(dots$adj, 2)
			ax[] <- adj[1]
			ay[] <- adj[2]
		}

		xmin <- xx - ax * w
		xmax <- xx + (1 - ax) * w
		ymin <- yy - ay * h
		ymax <- yy + (1 - ay) * h
		cx <- xy[, 1]
		cy <- xy[, 2]
		dx <- cbind(xmin - cx, xmax - cx, xmax - cx, xmin - cx)
		dy <- cbind(ymin - cy, ymin - cy, ymax - cy, ymax - cy)
		cs <- cos(srt * pi / 180)
		sn <- sin(srt * pi / 180)
		rx <- cx + dx * cs - dy * sn
		ry <- cy + dx * sn + dy * cs
		rx <- cbind(rx, rx[, 1])
		ry <- cbind(ry, ry[, 1])
		g <- cbind(
			geom = rep(seq_len(n), each = 5),
			part = 1,
			x = as.vector(t(rx)),
			y = as.vector(t(ry)),
			hole = 0
		)
		vect(g, type="polygons", crs="local")
	}


	.label_inside <- function(x, xy, labels, dots, halo=FALSE, hw=0.1) {
		n <- nrow(x)
		if (n == 0) {
			return(logical(0))
		}
		b <- try(.label_boxes(xy, labels, dots, halo, hw), silent=TRUE)
		if (inherits(b, "try-error") || nrow(b) != n) {
			return(rep(FALSE, n))
		}
		xx <- x
		crs(xx) <- "local"
		vapply(seq_len(n), function(i) {
			isTRUE(is.related(xx[i], b[i], "contains"))
		}, logical(1))
	}


	.central_part <- function(b) {
		if (nrow(b) == 0) {
			return(b)
		}
		b <- disagg(b)
		if (nrow(b) > 1) {
			b <- b[which.max(expanse(b, transform=FALSE))]
		}
		b
	}

	n <- nrow(x)
	if (n == 0) {
		return(cbind(x=numeric(), y=numeric()))
	}
	
	# same starting point as text(..., fit=FALSE)
	xy <- crds(centroids(x, inside=inside))
	# only shrink-search polygons whose drawn label does not fit
	if (grDevices::dev.cur() > 1) {
		need <- !.label_inside(x, xy, labels, dots, halo=halo, hw=hw)
	} else {
		need <- rep(TRUE, n)
	}
	if (!any(need)) {
		return(xy)
	}
	allxy <- xy
	x <- x[need]
	xy <- xy[need, , drop=FALSE]
	n <- nrow(x)
	# geodesic buffer is slow; treat coordinates as planar
	crs(x) <- "local"
	rad <- sqrt(expanse(x, transform=FALSE) / pi)
	# dense outlines: one light simplify makes every later buffer cheaper
	nv <- try(nrow(geom(x)), silent=TRUE)
	if (!inherits(nv, "try-error") && nv > 400) {
		tol <- stats::median(rad[is.finite(rad) & rad > 0], na.rm=TRUE) * 0.01
		if (is.finite(tol) && tol > 0) {
			xs <- try(simplifyGeom(x, tol), silent=TRUE)
			if (!inherits(xs, "try-error") && nrow(xs) == n) {
				x <- xs
				rad <- sqrt(expanse(x, transform=FALSE) / pi)
			}
		}
	}

	for (i in seq_len(n)) {
		p <- x[i]
		rmax <- rad[i]
		if (!is.finite(rmax) || rmax <= 0) {
			next
		}
		# holes make convHull(p2) ⊂ p impossible; test against the outer ring
		pouter <- fillHoles(p)
		# already convex (enough): true centroid if it lies inside
		if (is.related(pouter, hull(p, "convex"), "contains")) {
			cand <- centroids(p)
			if (is.related(p, cand, "contains")) {
				xy[i, ] <- crds(cand)
				next
			}
		}
		# binary search: smallest r whose inner hull fits, else the
		# largest non-empty shrink. Label placement only needs ~5% of rad.
		lo <- 0
		hi <- rmax
		tol <- rmax * 0.05
		found <- NULL
		last <- NULL
		for (k in 1:8) {
			if ((hi - lo) <= tol) {
				break
			}
			mid <- (lo + hi) / 2
			b <- buffer(p, -mid, quadsegs=3)
			a <- if (nrow(b) == 0) 0 else expanse(b, transform=FALSE)
			if (length(a) == 0 || is.na(a[1]) || a[1] <= 0) {
				hi <- mid
				next
			}
			b <- .central_part(b)
			last <- b
			if (is.related(pouter, hull(b, "convex"), "contains")) {
				found <- b
				hi <- mid
			} else {
				lo <- mid
			}
		}
		inner <- if (!is.null(found)) found else last
		if (!is.null(inner) && nrow(inner) > 0) {
			cand <- centroids(inner)
			if (!is.related(p, cand, "contains")) {
				cand <- centroids(inner, inside=TRUE)
			}
			if (is.related(p, cand, "contains")) {
				xy[i, ] <- crds(cand)
			}
		}
		if (!is.related(p, vect(xy[i, , drop=FALSE]), "contains")) {
			xy[i, ] <- crds(centroids(p, inside=TRUE))
		}
	}
	allxy[need, ] <- xy
	allxy
}



setMethod("text", signature(x="SpatVector"),
	function(x, labels, halo=FALSE, inside=FALSE, fit=FALSE, hc="white", hw=0.1, jitter=0, xpd=TRUE, ...) {
		if (missing(labels)) {
			labels <- 1:nrow(x)
		} else if (length(labels) == 1) {
			if (nrow(x) > 1) {
				labels <- as.data.frame(x)[,labels]
			} else {
				if (is.numeric(labels)) {
					if (labels %in% 1:ncol(x)) {
						labels <- x[[labels]][,1]
					}
				} else if (labels %in% names(x)) {
					labels <- x[[labels]][,1]
				}
			}
		}
		labels <- as.character(labels)
		if (length(labels) < nrow(x)) {
			labels <- rep_len(labels, nrow(x))
		}
		if (fit && (geomtype(x) == "polygons")) {
			xy <- .central_loc(x, labels, list(...), halo=halo, hw=hw, inside=inside)
		} else {
			xy <- crds(centroids(x, inside=inside))
		}
        if (isTRUE(jitter > 0)) {
			xy <- jitter(xy, factor = jitter)
		}
		if (halo && (isTRUE(hw > 0))) {
			.halo(xy[,1], xy[,2], labels, hc=hc, hw=hw, xpd=xpd, ...)
		} else {
			text(xy[,1], xy[,2], labels, xpd=xpd, ...)
		}
		invisible(xy)
	}
)



setMethod("boxplot", signature(x="SpatRaster"),
	function(x, y=NULL, maxcell=100000, ...) {
		if (is.null(y)) {
			cn <- names(x)
			if ( ncell(x) > maxcell) {
				warn("boxplot", "taking a sample of ", maxcell, " cells")
				x <- spatSample(x, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
			}
			names(x) <- cn
			boxplot(values(x), ...)
		} else {
			s <- c(x[[1]], y[[1]])
			if ( ncell(x) > maxcell) {
				warn("boxplot", "taking a regular sample of ", maxcell, " cells")
				s <- spatSample(s, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
			}
			s <- values(s, dataframe=TRUE)
			cn <- colnames(s)
			if (is.null(cn)) cn <- c("", "")
			colnames(s)[cn==""] <- c("layer1", "layer2")[cn==""]
			f <- try(stats::as.formula(paste(cn[1], '~', cn[2])), silent=TRUE)
			if (inherits(f, "try-error")) {
				colnames(s) <- c("layer1", "layer2")
				f <- layer1 ~ layer2
			}
			boxplot(f, data=s, ...)
		}
	}
)



setMethod("barplot", "SpatRaster",
	function(height, maxcell=1000000, digits=0, breaks=NULL, col, ...) {

		if (missing(col)) {
			col=grDevices::rainbow
		}
		height <- height[[1]]
		f <- is.factor(height)
#		x <- spatSample(height[[1]], maxcell, method="regular", as.raster=FALSE, as.df=f)
		x <- spatSample(height[[1]], maxcell, method="regular", as.raster=FALSE, as.df=FALSE, warn=FALSE)
		adj <- nrow(x) / ncell(height)
		if (adj < 1) {
			warn("barplot", "a sample of ", round(100*adj, 1), "% of the raster cells were used to estimate frequencies")
		}
		if (!f) {
			if (!is.null(digits)) {
				x <- round(x, digits)
			}
			if (!is.null(breaks)) {
				x <- cut(x, breaks)
			}
		}
		x <- table(x) / adj
		if (is.function(col)) {
			col <- col(length(x))
		}
		barplot(x, col=col, ...)
	}
)



shade <- function(slope, aspect, angle=45, direction=0, normalize=FALSE, filename="", overwrite=FALSE, ...) {
	stopifnot(inherits(slope, "SpatRaster"))
	opt <- spatOptions(filename, overwrite=overwrite, ...)
	slope@pntr <- slope@pntr$hillshade(aspect@pntr, angle, direction, normalize[1], opt)
	messages(slope, "shade")
}


map.pal <- function(name, n=50, ...) { 
	n <- round(n)
	if (n < 1) {
		error("map.pal", "n should be >= 1")	
	}
	f <- system.file("colors/palettes.rds", package="terra")
	v <- readRDS(f)
	if (missing(name)) {
		return(names(v))
	}
	if (name %in% names(v)) {
		v <- v[[name]]
		if ((n > 0) && ((length(v) != n))) {
			grDevices::colorRampPalette(v, ...)(n)
		} else {
			v
		}
	} else if (name == "random") {
		if (n > 433) {
			warning("map.pal", "you cannot get > 433 random colors; using n=433 instead")
			n <- 433
		}
		s <- sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)], n)
		rgb(t(grDevices::col2rgb(s))/255)

	} else {
		error("map.pal", paste(name, "is not a known palette"))
	}
}

map.leg <- function(name) { 
	f <- system.file("colors/legends.rds", package="terra")
	v <- readRDS(f)
	if (name %in% names(v)) {
		v[[name]]
	} else {
		error("map.leg", paste(name, "is not a known legend"))
	}
}

