# Author: Robert J. Hijmans
# Date :  October 2017
# Version 1.0
# License GPL v3


new_rast <- function(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=1, ymin=0, ymax=1, crs, extent, resolution, vals, names, time, units) {
		if (missing(extent)) {
			e <- c(xmin, xmax, ymin, ymax) 
		} else {
			e <- as.vector(extent)
		}
		if ((e[1] >= e[2]) || e[3] >= e[4]) {
			error("rast,missing", "invalid extent")
		}

		if (missing(crs)) {
			if (e[1] > -360.01 & e[2] < 360.01 & e[3] > -90.01 & e[4] < 90.01) {
				crs <- "OGC:CRS84"
			} else {
				crs <- ""
			}
		} else {
			crs <- as.character(crs)
		}
		#check_proj4_datum(crs)

		r <- methods::new("SpatRaster")
		r@ptr <- SpatRaster$new(c(nrows, ncols, nlyrs), e, crs)
		r <- messages(r, "rast")

		if (!missing(resolution)) {
			res(r) <- resolution
		}
		if (!missing(names)) {
			names(r) <- names
		}
		if (!missing(vals)) {
			values(r) <- vals
		}
		if (!missing(time)) {
			time(r) <- time
		}
		if (!missing(units)) {
			time(r) <- units
		}
		r
}


setMethod("rast", signature(x="missing"),
	function(x, nrows=180, ncols=360, nlyrs=1, xmin=-180, xmax=180, ymin=-90, ymax=90, crs, extent, resolution, vals, names, time, units) {
		new_rast(nrows, ncols, nlyrs, xmin, xmax, ymin, ymax, crs, extent, resolution, vals, names, time, units)
	}
)


setMethod("rast", signature(x="stars"),
	function(x) {
		x <- from_stars(x)
		if (inherits(x, "SpatRasterDataset")) {
			rast(x)
		} else {
			x
		}
	}
)

setMethod("rast", signature(x="stars_proxy"),
	function(x) {
		x <- from_stars(x)
		if (inherits(x, "SpatRasterDataset")) {
			rast(x)
		} else {
			x
		}
	}
)

setMethod("rast", signature(x="list"),
	function(x) {
		i <- sapply(x, function(i) inherits(i, "SpatRaster"))
		if (!all(i)) {
			if (!any(i)) {
				error("rast,list", "none of the elements of x are a SpatRaster")
			} else {
				warn("rast", sum(!i), " out of ", length(x), " elements of x are not a SpatRaster")
				x <- x[i]
			}
		}
		# start with an empty raster (alternatively use a deep copy)
		out <- rast(x[[1]])
		for (i in 1:length(x)) {
			out@ptr$addSource(x[[i]]@ptr, FALSE)
		}
		out <- messages(out, "rast")
		lnms <- names(x)
		i <- lnms != ""
		if (any(i)) {
			rnms <- names(out)
			rnms[lnms != ""] <- lnms[lnms != ""]
			names(out) <- rnms
		}
		out
	}
)



setMethod("rast", signature(x="SpatExtent"),
	function(x, ...) {
		dots <- list(...)
		dots$xmin=x[1]
		dots$xmax=x[2]
		dots$ymin=x[3]
		dots$ymax=x[4]
		do.call(new_rast, dots)
	}
)


setMethod("rast", signature(x="SpatVector"),
	function(x, ...) {
		dots <- list(...)
		e <- ext(x)
		dots$xmin=e[1]
		dots$xmax=e[2]
		dots$ymin=e[3]
		dots$ymax=e[4]
		if (all(is.na(pmatch(names(dots), "crs")))) {
			dots$crs <- crs(x)
		}
		do.call(new_rast, dots)
	}
)



.fullFilename <- function(x, mustExist=TRUE) {
	x <- trimws(x)
	p <- normalizePath(x, winslash = "/", mustWork = FALSE)
	if (mustExist) {
		i <- file.exists(p)
		if (all(i)) {
			return(p)
		} else {
			x[i] <- p[i]
		}
	} else {
		return(p)
	}
	#if (identical(basename(x), x)) {
	#	x <- file.path(getwd(), x)
	#}
	#if (expand) {
	#	x <- path.expand(x)
	#}
	return(x)
}

setMethod("rast", signature(x="character"),
	function(x, subds=0, lyrs=NULL, opts=NULL) {

		x <- trimws(x)
		x <- x[x!=""]
		if (length(x) == 0) {
			error("rast", "filename is empty. Provide a valid filename")
		}
		r <- methods::new("SpatRaster")
		f <- .fullFilename(x)
		#subds <- subds[1]
		if (is.null(opts)) opts <- ""[0]
		if (length(subds) == 0) subds = 0
		if (is.character(subds)) { 
			#r@ptr <- SpatRaster$new(f, -1, subds, FALSE, 0[])
			r@ptr <- SpatRaster$new(f, -1, subds, FALSE, opts, 0[])
		} else {
			r@ptr <- SpatRaster$new(f, subds-1, "", FALSE, opts, 0[])
		}
		r <- messages(r, "rast")
		if (r@ptr$getMessage() == "ncdf extent") {
			test <- try(r <- .ncdf_extent(r), silent=TRUE)
			if (inherits(test, "try-error")) {
				warn("rast", "GDAL did not find an extent. Cells not equally spaced?") 
			}
		}
		r <- messages(r, "rast")

		if (crs(r) == "") {
			if (is.lonlat(r, perhaps=TRUE, warn=FALSE)) {
				crs(r) <- "OGC:CRS84"
			}
		}
		
		if (!is.null(lyrs)) {
			r[[lyrs]]
		} else {
			r
		}
		
	}
)


multi <- function(x, subds=0, xyz=c(1,2,3)) {

	x <- trimws(x)
	x <- x[x!=""]
	if (length(x) == 0) {
		error("rast,character", "provide a valid filename")
	}
	r <- methods::new("SpatRaster")
	f <- .fullFilename(x)
	#subds <- subds[1]
	if (is.character(subds)) { 
		r@ptr <- SpatRaster$new(f, -1, subds, TRUE, xyz-1)
	} else {
		r@ptr <- SpatRaster$new(f, subds-1, "", TRUE, xyz-1)
	}
	if (r@ptr$getMessage() == "ncdf extent") {
		test <- try(r <- .ncdf_extent(r), silent=TRUE)
		if (inherits(test, "try-error")) {
			warn("rast", "GDAL did not find an extent. Cells not equally spaced?") 
		}
	}
	r <- messages(r, "rast")

	if (crs(r) == "") {
		if (is.lonlat(r, perhaps=TRUE, warn=FALSE)) {
			crs(r) <- "OGC:CRS84"
		}
	}
	r
}


setMethod("rast", signature(x="SpatRaster"),
	function(x, nlyrs=nlyr(x), names, vals, keeptime=TRUE, keepunits=FALSE, props=FALSE) {
		x@ptr <- x@ptr$geometry(nlyrs, props, keeptime, keepunits)
		x <- messages(x, "rast")
		if (!missing(names)) {
			if (length(names) == nlyr(x)) names(x) <- names
		}
		if (!missing(vals)) {
			values(x) <- vals
		}
		x
	}
)


setMethod("rast", signature(x="SpatRasterDataset"),
	function(x) {
		if (length(x) == 0) {
			error("rast", "empty SpatRasterDataset")
		} else if (length(x) == 1) {
			x[1]
		} else {
			r <- methods::new("SpatRaster")
			r@ptr <- x@ptr$collapse()
			nms <- names(x)
			if (any(nms != "")) {
				names(r) <- paste(rep(nms, nlyr(x)), names(r), sep="_")
			}
			r
		}
	}
)


setMethod("rast", signature(x="array"),
	function(x, ...) {
		dims <- dim(x)
		if (length(dims) > 3) {
			error("rast,array", "cannot handle an array with more than 3 dimensions")
		}
		r <- methods::new("SpatRaster")
		r@ptr <- SpatRaster$new(dims, c(0, dims[2], 0, dims[1]), "")
		values(r) <- x
		messages(r, "rast")
	}
)



setMethod("rast", signature(x="ANY"),
	function(x, ...) {
		methods::as(x, "SpatRaster")
	}
)


.rastFromXYZ <- function(xyz, digits=6, crs="", extent=NULL) {

	if (!is.null(e)) {
		warn("rast", 'argument "extent" is ignored if type="xyz"')
	}

	ln <- colnames(xyz)
	## xyz might not have colnames, or might have "" names
	if (is.null(ln)) ln <- rep("", ncol(xyz))
	if (any(nchar(ln) < 1)) ln <- make.names(ln)
	if (inherits(xyz, "data.frame")) {
		xyz <- as.matrix(xyz)
		xyz <- matrix(as.numeric(xyz), ncol=ncol(xyz), nrow=nrow(xyz))
	}
	x <- sort(unique(xyz[,1]))
	dx <- x[-1] - x[-length(x)]

	rx <- min(dx)
	for (i in 1:5) {
		rx <- rx / i
		q <- sum(round(dx / rx, digits=digits) %% 1)
		if ( q == 0 ) {
			break
		}
	}
	if ( q > 0 ) {
		error("raster,matrix(xyz)", "x cell sizes are not regular")
	}

	y <- sort(unique(xyz[,2]))
	dy <- y[-1] - y[-length(y)]
	# probably a mistake to use the line below
	# Gareth Davies suggested that it be removed
	# dy <- round(dy, digits)

	ry <- min(dy)
	for (i in 1:5) {
		ry <- ry / i
		q <- sum(round(dy / ry, digits=digits) %% 1)
		if ( q == 0 ) {
			break
		}
	}
	if ( q > 0 ) {
		error("raster,matrix(xyz)", "y cell sizes are not regular")
	}

	minx <- min(x) - 0.5 * rx
	maxx <- max(x) + 0.5 * rx
	miny <- min(y) - 0.5 * ry
	maxy <- max(y) + 0.5 * ry

	d <- dim(xyz)
	r <- rast(xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, crs=crs, nlyrs=d[2]-2)
	res(r) <- c(rx, ry)
	ext(r) <- round(ext(r), digits+2)
	cells <- cellFromXY(r, xyz[,1:2])
	if (d[2] > 2) {
		names(r) <- ln[-c(1:2)]
		v <- matrix(NA, nrow=ncell(r), ncol= nlyr(r))
		v[cells, ] <- xyz[, -c(1:2)]
		values(r) <- v
	}
	return(r)
}




setMethod("rast", signature(x="matrix"),
	function(x, type="", crs="", digits=6, extent=NULL) {
		stopifnot(prod(dim(x)) > 0)
		if (type == "xyz") {
			r <- .rastFromXYZ(x, crs=crs, digits=digits, extent=extent)
		} else {
			if (is.null(extent)) {
				r <- rast(nrows=nrow(x), ncols=ncol(x), crs=crs, extent=ext(c(0, 1, 0, 1)))
			} else {
				r <- rast(nrows=nrow(x), ncols=ncol(x), crs=crs, extent=extent)			
			}
			values(r) <- as.vector(t(x))
		}
		messages(r, "rast")
	}
)


setMethod("rast", signature(x="data.frame"),
	function(x, type="", crs="", digits=6, extent=NULL) {
		if (type == "xyz") {
			r <- .rastFromXYZ(x, crs=crs, digits=digits, extent=extent)
		} else {
			rast(as.matrix(x), type=type, crs=crs, digits=digits, extent=extent)
		}
	}
)


setMethod("NAflag<-", signature(x="SpatRaster"), 
	function(x, value)  {
		value <- as.numeric(value)
		if (!(x@ptr$setNAflag(value))) {
			error("NAflag<-", "cannot set this value")
		}
		x
	}
)

setMethod("NAflag", signature(x="SpatRaster"), 
	function(x)  {
		x@ptr$getNAflag()
	}
)

