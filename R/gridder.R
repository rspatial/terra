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
		if (any(!sapply(z, is.numeric))) {
			error(caller, paste("fields must be numeric"))
		}
	} else {
		if (!is.numeric(field)) {
			error(caller, paste(field, "is not numeric"))
		}
		z <- rep_len(field, nrow(y))
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


rastWinR <- function(x, y, win, pars, rfun, nl, filename, ...) {

	out <- rast(x, nlyr=nl)
	rb <- rast(out)
	e <- ext(out)
	hy <- yres(out)/2

	opt <- spatOptions()
	b <- writeStart(out, filename, n=12, sources="", ...)
	if (ncol(y) == 3) {
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
				p <- aggregate(p[2], p[1], rfun)
				v <- matrix(pars[5], ncell(rbe), nl)
				v[p[,1], ] <- p[,-1]
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
				p <- aggregate(y[p[[2]],3:ncol(y)], p[1], rfun)
				v <- matrix(pars[5], ncell(rbe), nl)
				v[p[,1], ] <- as.matrix(p[,-1])
				writeValues(out, v, b$row[i], b$nrows[i])					
			} else {
				writeValues(out, rep(pars[5], ncell(rbe) * nl), b$row[i], b$nrows[i])
			}
		}
	}
	writeStop(out)
}


rastBufR <- function(x, y, win, pars, rfun, nl, filename, ...) {
	w <- pars[1]
	z <- y[,3]
	y <- vect(y[,1:2,drop=FALSE])
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
	out <- rast(x, nlyr=1)
	ncs <- ncol(out)
	e <- ext(out)
	hy <- yres(out)/2
	off <- 0
	b <- writeStart(out, filename, n=12, sources="", ...)
	for (i in 1:b$n) {
		e$ymax <- yFromRow(out, b$row[i]) + hy
		e$ymin <- yFromRow(out, b$row[i] + b$nrows[i] - 1) - hy
		rbe <- crop(rb, e)
		f <- as.polygons(rbe, dissolve=FALSE)
		if (nrow(f) > 0) {
			f <- buffer(f, w)
			r <- relate(f, v, "intersects", pairs=TRUE)
		} else {
			r <- cbind(0,0)[0,]
		}
		if ((pars[4] > 1) && (nrow(r) > 0)) {
			a <- aggregate(r[,1], list(r[,1]), length)
			a <- a[a[,2] >= pars[4], 1]
			r <- r[r[,1] %in% a, ]
		} 
		if (nrow(r) > 0) {
			f <- aggregate(z[r[,2]], list(r[,1]), rfun)
			v <- matrix(pars[5], ncell(rbe), nl)
			v[f[,1], ] <- f[,-1]
			writeValues(out, v, b$row[i], b$nrows[i])					
		} else if (!is.na(pars[5])) {
			writeValues(out, rep(pars[5], ncell(rbe) * nl), b$row[i], b$nrows[i])
		}
	}
	writeStop(out)
}


setMethod("rasterizeWin", signature(x="SpatRaster", y="matrix"),
	function(x, y, fun, win="circle", pars, minPoints=1, fill=NA, filename="", ...) {

		pars <- c(get_rad(pars), minPoints[1], fill[1])
		if (ncol(y) < 3) {
			error("rasterizeNGB", "expecting a matrix with at least three columns")
		}
		win <- match.arg(tolower(win), c("circle", "ellipse", "rectangle", "buffer"))

		if (inherits(fun, "character")) {
			if (fun[1] == "count") {
				fun <- length
				pars[5] = 0
			}
			nl <- 1
		} else {
			i <- min(10, nrow(y))
			aggregate(y[1:i, -c(1:2)], list(rep(1, i)), fun)
			
			test <- fun(1:5)
			nl <- length(test)
		}
		nl <- nl * (ncol(y)-2)
		
		if (win == "buffer") {
			rastBufR(x, y, win, pars=pars, rfun=fun, nl=nl, filename=filename, ...)
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
					opt <- spatOptions(filename, ...)
					x@ptr <- x@ptr$rasterizeWindow(y[,1], y[,2], y[,3], algo, pars, opt)
					return(messages(x, "rasterizeWin"))
				}
			} 
			rastWinR(x, y, win, pars=pars, rfun=fun, nl=nl, filename=filename, ...)
		}
	}
)

setMethod("rasterizeWin", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, fun, win="circle", pars, minPoints=1, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("rasterizeWin", "SpatVector y must have a point geometry")		
		}
		y <- cbind(crds(y), get_z(y, field, "rasterizeWin"))
		rasterizeWin(x, y, fun=fun, win=win, pars=pars, minPoints=minPoints, fill=fill, filename=filename, ...)
	}
)




setMethod("interpNear", signature(x="SpatRaster", y="matrix"),
	function(x, y, radius, interpolate=FALSE, fill=NA, filename="", ...) {

		if (ncol(y) != 3) {
			error("interpNear", "expecting a matrix with three columns")
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
		interpNear(x, y, radius=radius, interpolate=interpolate, fill=fill, filename=filename, ...)
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="matrix"),
	function(x, y, radius, power=2, smooth=0, maxPoints=Inf, minPoints=1, near=FALSE, fill=NA, filename="", ...) {

		if (near) {
			algo <- "invdistpownear"
			pars <- c(power, smooth, radius[1], maxPoints, minPoints, fill)		
		} else {
			algo <- "invdistpow"
			pars <- c(power, smooth, get_rad(radius, "interpIDW"), maxPoints, minPoints, fill)			
		}

		if (ncol(y) != 3) {
			error("interpIDW", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rasterizeWindow(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpIDW")
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, radius, power=2, smooth=0, maxPoints=Inf, minPoints=1, near=FALSE, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpIDW", "SpatVector y must have a point geometry")		
		}
		y <- cbind(crds(y), get_z(y, field, "interpIDW"))
		interpIDW(x, y, radius, power=power, smooth=smooth, maxPoints=maxPoints, minPoints=minPoints, near=near, fill=fill, filename=filename, ...)
	}
)

