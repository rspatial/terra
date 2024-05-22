


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



.halo <- function(x, y=NULL, labels, col="black", hc="white", hw=0.1, ... ) {
	xy <- grDevices::xy.coords(x, y)
	xo <- hw * graphics::strwidth("A")
	yo <- hw * graphics::strheight("A")
	n <- nchar(labels)
	theta <- seq(pi/4, 2*pi, length.out=8*hw*10)
	for (i in theta) {
		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=hc, ... )
	}
	graphics::text(xy$x, xy$y, labels, col=col, ... )
}

halo <- function(x, y=NULL, labels, col="black", hc="white", hw=0.1, ... ) {
	.halo(x=x, y=y, labels=labels, col=col, hc=hc, hw=hw, ...) 
}


setMethod("text", signature(x="SpatRaster"),
	function(x, labels, digits=0, halo=FALSE, ...) {
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
		if (is.factor(labels)) {
			labels <- substr(as.character(labels), 1, max(1, digits))
		} else if (is.numeric(labels)) {
			labels <- as.character(round(labels, digits=digits) )
		}
		if (halo) {
			.halo(xy[,1], xy[,2], labels, ...)
		} else {
			text(xy[,1], xy[,2], labels, ...)
		}
	}
)


setMethod("text", signature(x="SpatVector"),
	function(x, labels, halo=FALSE, inside=FALSE, ...) {
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
		xy <- geom(centroids(x, inside=inside))[,c("x","y"),drop=FALSE]
		if (halo) {
			.halo(xy[,1], xy[,2], labels, ...)
		} else {
			text(xy[,1], xy[,2], labels, ...)
		}
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
	slope@ptr <- slope@ptr$hillshade(aspect@ptr, angle, direction, normalize[1], opt)
	messages(slope, "shade")
}


map.pal <- function(name, n=50, ...) { 
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

