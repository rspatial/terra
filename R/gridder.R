
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
		if (length(field) != 1) {
			error(caller, "field name should have length 1")
		}
		if (!field %in% names(y)) {
			error(caller, paste(field, "is not a name in y"))			
		}
		z <- y[[field, drop=TRUE]]
		if (!is.numeric(z)) {
			error(caller, paste(field, "is not numeric"))
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
		pars <- c(get_rad(radius), minPoints[1], fill[1])
		if (ncol(y) != 3) {
			error("rasterizeNGB", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		if (builtin) {
			x@ptr <- x@ptr$rasterizeWindow(y[,1], y[,2], y[,3], algo, pars, opt)
			messages(x, "rasterizeWin")
		} else {
			p <- x@ptr$winpoints(y[,1], y[,2], pars, opt)
			x <- rast(x, nlyr=1)
			if (!is.na(fill[1])) {
				x <- init(x, fill[1])
			}
			
			if (length(p[[1]]) == 0) {
				warn("rasterizeWin", "All windows were empty")
			} else {
				p <- aggregate(list(y[,3][p[[2]]+1]), p[1], fun)
				set.values(x, p[,1]+1, p[,2])
			}
			x
		}
	}
)

setMethod("rasterizeWin", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, fun, radius, minPoints=1, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("rasterizeWin", "SpatVector y must have a point geometry")		
		}
		y <- cbind(crds(y), get_z(y, field, "rasterizeWin"))
		rasterizeWin(x, y, fun=fun, radius=radius, minPoints=minPoints, fill=fill, filename=filename, ...)
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

