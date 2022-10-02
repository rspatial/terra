
check_ngb_pars <- function(algo, pars, fill) {
	#p <- c("Power", "Smoothing", "Radius", "Radius1", "Radius2", "Angle", "nMaxPoints", "nMinPoints")
	n <- tolower(names(pars))
	if (length(n) == 0) error("rasterizeNGB", "parameters are not named") 	
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
		error("rasterizeNGB", "invalid algorithm name")
	}
	if (length(n) != length(pex)) error("rasterizeNGB", paste("expected 4 parameters for", algo, "got",  length(n)))
	if (!all(n %in% tolower(pex)))error("rasterizeNGB", paste("parameters needed for", algo, ":\n  ", paste(pex, collapse=",")))
	c(pars[match(n, tolower(pex))], fill)
}


get_z <- function(y, field) {
	if (inherits(field, "character")) {
		if (length(field != 1)) {
			error("rasterizeNGB", "field name should have length 1")
		}
		if (!field %in% names(y)) {
			error("rasterizeNGB", paste(field, "is not a name in y"))			
		}
		z <- y[[field, drop=TRUE]]
		if (!is.numeric(z)) {
			error("rasterizeNGB", paste(field, "is not numeric"))
		}
	} else {
		if (!is.numeric(field)) {
			error("rasterizeNGB", paste(field, "is not numeric"))
		}
		z <- rep_len(field, nrow(y))
	}
	z
}

get_rad <- function(r, caller="rasterizeNGB") {
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



setMethod("rasterizeWin", signature(x="SpatRaster", y="matrix"),
	function(x, y, fun, radius, minPoints=1, fill=NA, filename="", ...) {

		algo <- .makeTextFun(fun)
		algos <- c("min", "max", "range", "mean", "count", "distto", "distbetween")
		builtin <- FALSE
		if (inherits(algo, "character")) {
			if ((algo %in% algos)) {
				builtin <- TRUE
			}
		} 
		pars <- c(get_rad(radius), minPoints, fill)
		if (ncol(y) != 3) {
			error("rasterizeNGB", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		if (builtin) {
			x@ptr <- x@ptr$rasterizeWindow(y[,1], y[,2], y[,3], algo, pars, opt)
			messages(x, "rasterizeNGB")
		} else {		
			p <- x@ptr$winpoints(y[,1], y[,2], pars, opt)
			p <- aggregate(list(y[,3][p[[2]]+1]), p[1], fun)
			x <- rast(x, nlyr=1)
			x[p[,1]+1] <- p[,2]
			x
		}
	}
)

setMethod("rasterizeWin", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, fun, radius, minPoints=1, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("rasterizeNGB", "SpatVector y must have a point geometry")		
		}
		y <- cbind(crds(y), get_z(y, field))
		raterizeNGB(x, y, fun=fun, radius=radius, minPoints=minPoints, fill=fill, filename=filename, ...)
	}
)




setMethod("interpNear", signature(x="SpatRaster", y="matrix"),
	function(x, y, algo, radius, fill=NA, filename="", ...) {

		algos <- c("nearest", "linear")
		algo <- match.arg(tolower(algo), algos)

		if (algo == "nearest") {
			pars <- c(get_rad(radius, "interpNear"), fill)
		} else {
			pars <- c(radius[1], fill)		
		}

		if (ncol(y) != 3) {
			error("interpNear", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rasterizeWin(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpNear")
	}
)

setMethod("interpNear", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, algo, radius, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpNear", "SpatVector y must have a point geometry")		
		}
		y <- cbind(crds(y), get_z(y, field))
		interpNear(x, y, algo=algo, radius=radius, fill=fill, filename=filename, ...)
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="matrix"),
	function(x, y, algo, radius, power=2, smooth=0, maxPoints, minPoints=1, near=FALSE, fill=NA, filename="", ...) {

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
		x@ptr <- x@ptr$rasterizeWin(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpIDW")
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, radius, power=2, smooth=0, maxPoints=Inf, minPoints=1, near=FALSE, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpIDW", "SpatVector y must have a point geometry")		
		}
		y <- cbind(crds(y), get_z(y, field))
		interpIDW(x, y, radius, power=power, smooth=smooth, maxPoints=maxPoints, minPoints=minPoints, near=near, fill=fill, filename=filename, ...)
	}
)

