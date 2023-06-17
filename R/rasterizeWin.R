# Author: Robert J. Hijmans
# Date : October 2022
# Version 1.0
# License GPL v3


check_ngb_pars <- function(algo, pars, fill, caller="rasterizeWin") {
	#p <- c("Power", "Smoothing", "Radius", "Radius1", "Radius2", "Angle", "nMaxPoints", "nMinPoints")
	n <- tolower(names(pars))
	if (length(n) == 0) error(caller, "parameters are not named") 
	if (algo %in%  c("min", "max", "range", "mean", "count", "distto", "distbetween")) {
		pex <- c("Radius1", "Radius2", "Angle", "MinPoints")
	} else if (algo == "invdistpow") {
		pex <- c("Power", "Smoothing", "Radius1", "Radius2", "Angle", "MaxPoints", "MinPoints")
	} else if (algo == "invdistpownear") {
		pex <- c("Power", "Smoothing", "Radius", "MaxPoints", "MinPoints")
	} else if (algo == "nearest") {
		pex <- c("Radius1", "Radius2", "Angle")
	} else if (algo == "linear") {
		pex <- c("Radius")
	} else {
		error(caller, "invalid algorithm name")
	}
	if (length(n) != length(pex)) error(caller, paste("expected 4 parameters for", algo, "got",  length(n)))
	if (!all(n %in% tolower(pex)))error(caller, paste("parameters needed for", algo, ":\n  ", paste(pex, collapse=",")))
	c(pars[match(n, tolower(pex))], fill)
}


get_z <- function(y, field, caller="rasterizeWin") {
	if (inherits(field, "character")) {
		if (length(field) == 0) {
			error(caller, "field name missing")
		}
		if (!all(field %in% names(y))) {
			f <- paste(field[!(field %in% names(y))], collapse=", ")
			error(caller, paste(f, " is not a name in y"))
		}
		z <- y[, field, drop=TRUE]
		#if (any(!sapply(z, is.numeric))) {
		#	error(caller, paste("fields must be numeric"))
		#}
	} else {
		#if (!is.numeric(field)) {
		#	error(caller, paste(field, "is not numeric"))
		#}
		z <- data.frame(field=rep_len(field, nrow(y)))
	}
	z
}

get_rad <- function(r, caller="rasterizeWin") {
	if (length(r) == 1) {
		c(r, r, 0)
	} else if (length(r) == 2) {
		c(r, 0)
	} else if (length(r) == 3) {
		r[3] = r[3] %% 360
		if (r[3] < 0) r[3] = r[3] + 360;
		r
	} else {
		error(caller, "radius argument should have length 1, 2, or 3")
	}
}


rastWinR <- function(x, y, win, pars, fun, nl, cvars, filename, wopt, ...) {

	out <- rast(x, nlyr=nl)
	rb <- rast(out)
	e <- ext(out)
	hy <- yres(out)/2

	opt <- spatOptions()
	b <- writeStart(out, filename, n=12, sources="", wopt=wopt)
	if ((ncol(y) == 3) && (is.numeric(y[,3]))) {
		for (i in 1:b$n) {
			e$ymax <- yFromRow(out, b$row[i]) + hy
			e$ymin <- yFromRow(out, b$row[i]  + b$nrows[i] - 1) - hy
			rbe <- crop(rb, e)
			if (win == "rectangle") {
				p <- rbe@ptr$winrect(y[,1], y[,2], y[,3], pars, opt)
			} else {
				p <- rbe@ptr$wincircle(y[,1], y[,2], y[,3], pars, opt)
			}

			if ((pars[4] > 1) && (length(p[[1]]) > 0)) {
				a <- aggregate(p[[1]], p[1], length)
				a <- a[a[,2] >= pars[4], 1]
				a <- p[[1]] %in% a
				p[[1]] <- p[[1]][a]
				p[[2]] <- p[[2]][a]
			} 

			if (length(p[[1]]) > 0) {
				p <- aggregate(p[2], p[1], fun, ...)
				v <- matrix(pars[5], ncell(rbe), nl)
				v[p[,1]+1, ] <- p[,-1]
				writeValues(out, v, b$row[i], b$nrows[i])
			} else {
				writeValues(out, rep(pars[5], ncell(rbe) * nl), b$row[i], b$nrows[i])
			}
		}
	} else {
		id <- 1:nrow(y)
		for (i in 1:b$n) {
			e$ymax <- yFromRow(out, b$row[i]) + hy
			e$ymin <- yFromRow(out, b$row[i]  + b$nrows[i] - 1) - hy
			rbe <- crop(rb, e)
			if (win == "rectangle") {
				p <- rbe@ptr$winrect(y[,1], y[,2], id, pars, opt)
			} else {
				p <- rbe@ptr$wincircle(y[,1], y[,2], id, pars, opt)
			}

			if ((pars[4] > 1) && (length(p[[1]]) > 0)) {
				a <- aggregate(p[[1]], p[1], length)
				a <- a[a[,2] >= pars[4], 1]
				a <- p[[1]] %in% a
				p[[1]] <- p[[1]][a]
				p[[2]] <- p[[2]][a]
			} 
	
			if (length(p[[1]]) > 0) {
				py <- y[p[[2]],3:ncol(y)]
				v <- matrix(pars[5], ncell(rbe), nl)
				if (cvars) {
					#u <- unique(p[[1]])
					#if (usedots) {
					#	p <- sapply(u, function(i) fun(py[p[[1]]==i, ,drop=FALSE], ...))
					#} else {
					#	p <- sapply(u, function(i) fun(py[p[[1]]==i, ,drop=FALSE]))
					#}
					p <- split(py, p[[1]])
					p <- sapply(p, fun, ...)
					if (!is.null(dim(p))) {
						p <- t(as.matrix(p))
						u <- rownames(p)
					} else {
						u <- names(p)
					}
					u <- as.numeric(u)
					v[u+1, ] <- p
				} else {
					p <- aggregate(py, p[1], fun, ...)
					v[p[,1]+1, ] <- as.matrix(p[,-1])
				}
				writeValues(out, v, b$row[i], b$nrows[i])
			} else {
				writeValues(out, rep(pars[5], ncell(rbe) * nl), b$row[i], b$nrows[i])
			}
		}
	}
	writeStop(out)
}


rastBufR <- function(x, y, win, pars, fun, nl, cvars, filename, wopt, ...) {
	w <- pars[1]
	z <- y[,-c(1:2), drop=FALSE]
	y <- vect(y[,1:2,drop=FALSE], geom=c("x", "y"))
	out <- rast(x, nlyr=1)
	rb <- rast(out)
	if (!is.lonlat(x)) {
		ngb <- max(round(w/yres(x)), round(w/xres(x)))
		if (ngb <= 5) {
			m <- matrix(1, 2*round(w/yres(x))+1, 2*round(w/xres(x))+1)
			rb <- rasterize(y, rb)
			rb <- focal(rb, m, "sum", na.rm=TRUE)
		} 
	}
	out <- rast(x, nlyr=nl)
	ncs <- ncol(out)
	e <- ext(out)
	hy <- yres(out)/2
	off <- 0
	b <- writeStart(out, filename, n=12, sources="", wopt=wopt)

	for (i in 1:b$n) {
		e$ymax <- yFromRow(out, b$row[i]) + hy
		e$ymin <- yFromRow(out, b$row[i] + b$nrows[i] - 1) - hy
		rbe <- crop(rb, e)
		buf <- as.polygons(rbe, dissolve=FALSE)
		buf <- buffer(buf, w)
		r <- relate(buf, y, "intersects", pairs=TRUE)
		if ((pars[4] > 1) && (nrow(r) > 0)) {
			a <- aggregate(r[,1], list(r[,1]), length)
			a <- a[a[,2] >= pars[4], 1]
			r <- r[r[,1] %in% a, ]
		} 
		if (nrow(r) > 0) {
			v <- matrix(pars[5], ncell(rbe), nl)
			if ((ncol(z) == 1) || (!cvars)) {
				f <- aggregate(z[r[,2],,drop=FALSE], list(r[,1]), fun, ...)
				v[f[,1], ] <- f[,-1]
				writeValues(out, v, b$row[i], b$nrows[i])
			} else {
				#u <- unique(r[,1])
				#p <- z[r[,2], ,drop=FALSE]
				#if (usedots) {
				#	p <- sapply(u, function(i) fun(p[r[,1]==i, ,drop=FALSE], ...))
				#} else {
				#	p <- sapply(u, function(i) fun(p[r[,1]==i, ,drop=FALSE]))
				#}
				#if (!is.null(dim(p))) p <- t(p)
				#v[u, ] <- as.matrix(p)

				py <- z[r[,2], ,drop=FALSE]
				s <- split(py, r[,1])
				p <- sapply(s, fun, ...)
				if (!is.null(dim(p))) {
					p <- t(as.matrix(p))
					u <- as.numeric(rownames(p))
				} else {
					u <- as.numeric(names(p))
				}
				v[u, ] <- p
				writeValues(out, v, b$row[i], b$nrows[i])
			}
		} else {
			writeValues(out, rep(pars[5], ncell(rbe) * nl), b$row[i], b$nrows[i])
		}
	}
	writeStop(out)
}


setMethod("rasterizeWin", signature(x="data.frame", y="SpatRaster"),
	function(x, y, win="circle", pars, fun, ..., cvars=FALSE, minPoints=1, fill=NA, filename="", wopt=list()) {

		pars <- c(get_rad(pars), minPoints[1], fill[1])
		if (ncol(x) < 3) {
			error("rasterizeNGB", "expecting a matrix with at least three columns")
		}
		win <- match.arg(tolower(win), c("circle", "ellipse", "rectangle", "buffer"))

#		usedots <- length(list(...)) > 0
		if (ncol(x) == 3) cvars = FALSE
		if (inherits(fun, "character")) {
			if (fun[1] == "count") {
				fun <- length
				pars[5] = 0
			}
			nl <- ncol(x)-2
		} else {
			if (cvars) {
				i <- min(10, nrow(x))
				v <- x[1:i, -c(1:2)]
				test <- sapply(list(v), fun, ...)
				nl <- length(test)
			} else {
				test <- sapply(list(1:5), fun, ...)
				nl <- length(test) * (ncol(x)-2)
			}
		}

		if (win == "buffer") {
			rastBufR(y, x, win, pars=pars, fun=fun, nl=nl, cvars=cvars, filename=filename, wopt=wopt, ...)
		} else {
			if (win == "circle") {
				pars[2] = pars[1]
				pars[3] = 0;
			}
			algo <- .makeTextFun(fun)
			#algos <- c("min", "max", "range", "mean", "count", "distto", "distbetween")
			algos <- c("distto", "distbetween")
			builtin <- FALSE
			if (inherits(algo, "character") && (algo %in% algos)) {
				if (win == "rectangle") {
					error("rasterizeWin", paste(fun, "not yet available for 'win=rectangle'"))
				} else {
					opt <- spatOptions(filename, wopt=wopt)
					x@ptr <- x@ptr$rasterizeWindow(x[,1], x[,2], x[,3], algo, pars, opt)
					return(messages(x, "rasterizeWin"))
				}
			} 
			rastWinR(x=y, y=x, win=win, pars=pars, fun=fun, nl=nl, cvars=cvars, filename=filename, wopt=wopt, ...)
		}
	}
)

setMethod("rasterizeWin", signature(x="SpatVector", y="SpatRaster"),
	function(x, y, field, win="circle", pars, fun, ..., cvars=FALSE, minPoints=1, fill=NA, filename="", wopt=list()) {
		if (geomtype(x) != "points") {
			error("rasterizeWin", "SpatVector y must have a point geometry")
		}
		x <- cbind(crds(x), get_z(x, field, "rasterizeWin"))
		rasterizeWin(x, y, win=win, pars=pars, fun=fun, minPoints=minPoints, fill=fill, cvars=cvars, filename=filename, wopt=wopt, ...)
	}
)




setMethod("interpNear", signature(x="SpatRaster", y="matrix"),
	function(x, y, radius, interpolate=FALSE, fill=NA, filename="", ...) {

		if (ncol(y) != 3) {
			error("interpNear", "expecting a matrix with three columns")
		}
		if (!is.numeric(y)) {
			error("interpNear", "values must be numeric")
		}

		if (interpolate) {
			algo <- "linear"
			pars <- c(radius[1], fill)
		} else {
			algo <- "nearest"
			pars <- c(get_rad(radius, "interpNear"), fill)
		}

		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rasterizeWindow(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpNear")
	}
)


setMethod("interpNear", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, radius, interpolate=FALSE, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpNear", "SpatVector y must have a point geometry")
		}
		y <- cbind(crds(y), get_z(y, field, "interpNear"))
		y <- as.matrix(y)
		interpNear(x, y, radius=radius, interpolate=interpolate, fill=fill, filename=filename, ...)
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="matrix"),
	function(x, y, radius, power=2, smooth=0, maxPoints=Inf, minPoints=1, near=TRUE, fill=NA, filename="", ...) {

		if (ncol(y) != 3) {
			error("interpIDW", "expecting a matrix with three columns")
		}
		if (!is.numeric(y)) {
			error("interpIDW", "values must be numeric")
		}

		if (near) {
			algo <- "invdistpownear"
			pars <- c(power, smooth, radius[1], maxPoints, minPoints, fill)
		} else {
			algo <- "invdistpow"
			pars <- c(power, smooth, get_rad(radius, "interpIDW"), maxPoints, minPoints, fill)
		}

		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rasterizeWindow(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpIDW")
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, radius, power=2, smooth=0, maxPoints=Inf, minPoints=1, near=TRUE, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpIDW", "SpatVector y must have a point geometry")
		}
		y <- cbind(crds(y), get_z(y, field, "interpIDW"))
		y <- as.matrix(y)
		interpIDW(x, y, radius, power=power, smooth=smooth, maxPoints=maxPoints, minPoints=minPoints, near=near, fill=fill, filename=filename, ...)
	}
)

