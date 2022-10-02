
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

get_rad <- function(r) {
	if (length(r) == 1) {
		c(r, r, 0)
	} else if (length(r) == 2) {
		c(r, 0)
	} else if (length(r) == 3) {
		r
	} else {
		error("rasterizeNGB", "radius should have length 1, 2, or 3")
	}
}

setMethod("rasterizeNGB", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, algo, radius, minPoints=1, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("rasterizeNGB", "SpatVector y must have a point geometry")		
		}
		algos <- c("min", "max", "range", "mean", "count", "distto", "distbetween")
		algo <- match.arg(tolower(algo), algos)		
		pars <- c(get_rad(radius), minPoints, fill)
		z <- get_z(y, field)
		y <- crds(y)
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$gridder(y[,1], y[,2], z, algo, pars, opt)
		messages(x, "rasterizeNGB")
	}
)



setMethod("rasterizeNGB", signature(x="SpatRaster", y="matrix"),
	function(x, y, algo, radius, minPoints=1, fill=NA, filename="", ...) {

		algos <- c("min", "max", "range", "mean", "count", "distto", "distbetween")
		algo <- match.arg(tolower(algo), algos)
		pars <- c(get_rad(radius), minPoints, fill)
		if (ncol(y) != 3) {
			error("rasterizeNGB", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$gridder(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "rasterizeNGB")
	}
)


setMethod("interpNear", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, algo, radius, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpNear", "SpatVector y must have a point geometry")		
		}
		algos <- c("nearest", "linear")
		algo <- match.arg(tolower(algo), algos)
		if (algo == "nearest") {
			pars <- c(get_rad(radius), fill)
		} else {
			pars <- c(radius[1], fill)		
		}

		z <- get_z(y, field)
		y <- crds(y)
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$gridder(y[,1], y[,2], z, algo, pars, opt)
		messages(x, "interpNear")
	}
)



setMethod("interpNear", signature(x="SpatRaster", y="matrix"),
	function(x, y, algo, radius, fill=NA, filename="", ...) {

		algos <- c("nearest", "linear")
		algo <- match.arg(tolower(algo), algos)

		if (algo == "nearest") {
			pars <- c(get_rad(radius), fill)
		} else {
			pars <- c(radius[1], fill)		
		}

		if (ncol(y) != 3) {
			error("interpNear", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$gridder(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpNear")
	}
)


setMethod("interpIDW", signature(x="SpatRaster", y="SpatVector"),
	function(x, y, field, algo, radius, power=2, smooth=0, maxPoints, minPoints=1, fill=NA, filename="", ...) {
		if (geomtype(y) != "points") {
			error("interpIDW", "SpatVector y must have a point geometry")		
		}
		algos <- c("invdistpow", "invdistpownear")
		algo <- match.arg(tolower(algo), algos)

		if (algo == "invdistpow") {
			pars <- c(power, smooth, get_rad(radius), maxPoints, minPoints, fill)
		} else {
			pars <- c(power, smooth, radius[1], maxPoints, minPoints, fill)		
		}
		z <- get_z(y, field)
		y <- crds(y)
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$gridder(y[,1], y[,2], z, algo, pars, opt)
		messages(x, "interpIDW")
	}
)



setMethod("interpIDW", signature(x="SpatRaster", y="matrix"),
	function(x, y, algo, radius, power=2, smooth=0, maxPoints, minPoints=1, fill=NA, filename="", ...) {

		algos <- c("invdistpow", "invdistpownear")
		algo <- match.arg(tolower(algo), algos)
		if (algo == "invdistpow") {
			pars <- c(power, smooth, get_rad(radius), maxPoints, minPoints, fill)
		} else {
			pars <- c(power, smooth, radius[1], maxPoints, minPoints, fill)		
		}

		if (ncol(y) != 3) {
			error("interpIDW", "expecting a matrix with three columns")
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$gridder(y[,1], y[,2], y[,3], algo, pars, opt)
		messages(x, "interpIDW")
	}
)
