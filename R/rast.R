# Author: Robert J. Hijmans
# Date :  October 2017
# Version 1.0
# License GPL v3


new_rast <- function(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=1, ymin=0, ymax=1, crs, extent, resolution, vals, names, time, units) {

	ncols <- round(ncols)
	if (ncols < 1) error("rast", "ncols < 1")
	nrows <- round(nrows)
	if (nrows < 1) error("rast", "nrows < 1")

	if (missing(extent)) {
		e <- c(xmin, xmax, ymin, ymax)
	} else {
		extent <- ext(extent)
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
		crs <- character_crs(crs, "rast")
	}
	#check_proj4_datum(crs)

	r <- methods::new("SpatRaster")
	r@ptr <- SpatRaster$new(c(nrows, ncols, nlyrs), e, crs)
	r <- messages(r, "rast")

	if (!missing(resolution)) {
		res(r) <- resolution
	}
	if (!missing(vals)) {
		if (length(vals) == 1) {
			if (is.na(vals[1])) {
				vals <- as.numeric(NA)
			}
		}
		if (!is.null(ncol(vals))) {
			nms <- colnames(vals)
			if (length(nms) == nlyr(r)) {
				names(r) <- nms
			}
		}
		values(r) <- vals
	}
	if (!missing(names)) {
		names(r) <- names
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




rast_from_image <- function(x) { 
	# list representing an "image"
	if (! all(dim(x$z) == c(length(x$x), length(x$y)))) {
		error("rast", '"z" does not have the right dimensions') 
	}
	# omitted "-1" bug fix by Barry Rowlingson 
	resx <- ( x$x[length(x$x)] - x$x[1] ) / (length(x$x)-1)
	resy <- ( x$y[length(x$y)] - x$y[1] ) / (length(x$y)-1)
	xmn <- min(x$x) - 0.5 * resx
	xmx <- max(x$x) + 0.5 * resx
	ymn <- min(x$y) - 0.5 * resy
	ymx <- max(x$y) + 0.5 * resy

	dx <- abs(max(abs((x$x[-1] - x$x[-length(x$x)])) / resx) - 1)
	dy <- abs(max(abs((x$y[-1] - x$y[-length(x$y)])) / resy) - 1)
	if (is.na(dx) | is.na(dy)) {
		error("rast", "missing values in coordinates")
	} 
	if (dx > 0.01 | dy > 0.01) {
		error("rast", "data are not on a regular grid")
	}
#	if (xmn > -360.1 & xmx < 360.1 & ymn > -90.1 & ymx < 90.1) { 
#		crs <- "+proj=longlat +datum=WGS84"
#	} else {
#		crs <- ""
#	}
	x <- t(x$z)
	x <- x[nrow(x):1, ]
	rast(x, extent=ext(xmn, xmx, ymn, ymx), crs="")
}



setMethod("rast", signature(x="list"),
	function(x, warn=TRUE) {
		i <- sapply(x, function(i) inherits(i, "SpatRaster"))
		if (!all(i)) {
			if (!any(i)) {
				if ((length(x) == 3) && all(c("x", "y", "z") %in% names(x))) {
					return(rast_from_image(x))
				}
				error("rast,list", "none of the elements of x are a SpatRaster")
			} else {
				warn("rast", sum(!i), " out of ", length(x), " elements of x are not a SpatRaster")
				x <- x[i]
			}
		}
		# start with an empty raster (alternatively use a deep copy)
		out <- deepcopy(x[[1]])
		if (length(x) == 1) {
			return(out)
		}
		opt <- spatOptions()
		for (i in 2:length(x)) {
			out@ptr$addSource(x[[i]]@ptr, warn, opt)
		}
		out <- messages(out, "rast")
		lnms <- names(x)
		if (!is.null(lnms)) {
			if (any(lnms != "") && (length(lnms) == nlyr(out))) {
				rnms <- names(out)
				rnms[lnms != ""] <- lnms[lnms != ""]
				names(out) <- rnms
			} else if (all(lnms != "")) {
				nl <- sapply(x, nlyr)
				rnms <- sapply(1:length(nl), function(i) {
							if (nl[i] > 1) paste0(lnms[i], "_", 1:nl[i]) else lnms[i]
						})
				names(out) <- unlist(rnms)
			}
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



.fullFilename <- function(x, mustExist=FALSE, vsi=FALSE) {
	x <- trimws(x)
	x <- x[x != ""]
	x <- enc2utf8(x)
	
	i <- substr(x, 1, 5) == "s3://" 
	x[i] <- paste0("/vsis3/", substr(x[i], 6, nchar(x[i])))
	if (all(i)) return(x)
	
	i <- substr(x, 1, 4) == "http" 
	if (vsi) {
		x[i] <- paste0("/vsicurl/", x[i])
	}
	if (all(i)) return(x)

	i <- grepl(":", x)
	if (all(i)) return(x)
	
	p <- normalizePath(x[!i], winslash = "/", mustWork = FALSE)
	if (mustExist) {
		j <- file.exists(dirname(p))
		x[j] <- p[j]
	} else {
		i <- !file.exists(p)
		p[i] <- x[i]
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
	function(x, subds=0, lyrs=NULL, drivers=NULL, opts=NULL, win=NULL, snap="near", vsi=FALSE, raw=FALSE) {

		f <- .fullFilename(x, vsi=vsi)
		if (length(f) == 0) {
			error("rast", "filename is empty. Provide a valid filename")
		}

		if ((length(f) == 1) && grepl("\\.rds$", tolower(f[1]))) {
			r <- unwrap(readRDS(x))
			if (!inherits(r, "SpatRaster")) {
				error("rast", "the rds file does not store a SpatRaster")
			}
			return(r)
		}
		
		r <- methods::new("SpatRaster")
		#subds <- subds[1]
		if (is.null(opts)) opts <- ""[0]
		if (raw) opts <- c(opts, "so=false")
		if (is.null(drivers)) drivers <- ""[0]
		if (length(subds) == 0) subds = 0
		if (is.character(subds)) {
			#r@ptr <- SpatRaster$new(f, -1, subds, FALSE, 0[])
			r@ptr <- SpatRaster$new(f, -1, subds, FALSE, drivers, opts, 0[])
		} else {
			r@ptr <- SpatRaster$new(f, subds-1, "", FALSE, drivers, opts, 0[])
		}
		r <- messages(r, "rast")
		if (r@ptr$getMessage() == "ncdf extent") {
			# could have used opts="IGNORE_XY_AXIS_NAME_CHECKS=YES"
			test <- try(r <- .ncdf_extent(r, f), silent=TRUE)
			if (inherits(test, "try-error")) {
				warn("rast", "GDAL did not find an extent. Cells not equally spaced?")
			}
		}
		r <- messages(r, "rast")
		if (crs(r) == "") {
			if (is.lonlat(r, perhaps=TRUE, warn=FALSE)) {
				if (!isTRUE(all(as.vector(ext(r)) == c(0,ncol(r),0,nrow(r))))) {
					crs(r) <- "OGC:CRS84"
				}
			}
		}

		if (!is.null(lyrs)) {
			r <- r[[lyrs]]
		} 
		if (!is.null(win)) {
			e <- ext(win)
			e <- align(e, r, snap=snap)
			window(r) <- e
		}
		r
	}
)


multi <- function(x, subds=0, xyz=3:1, drivers=NULL, opts=NULL) {

	x <- trimws(x)
	x <- x[x!=""]
	if (length(x) == 0) {
		error("rast,character", "provide a valid filename")
	}
	r <- methods::new("SpatRaster")
	f <- .fullFilename(x)
	if (is.null(opts)) opts <- ""[0]
	if (is.null(drivers)) drivers <- ""[0]
	if (length(subds) == 0) subds = 1
	subds <- subds[1]

	if (is.character(subds)) {
		r@ptr <- SpatRaster$new(f, -1, subds, TRUE, drivers, opts, xyz-1)
	} else {
		r@ptr <- SpatRaster$new(f, subds-1, ""[0], TRUE, drivers, opts, xyz-1)
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
	function(x, nlyrs=nlyr(x), names, vals, keeptime=TRUE, keepunits=FALSE, props=FALSE, tags=FALSE) {
		if (inherits(nlyrs, "SpatRaster")) {
			error("rast", "use 'c()' to combine SpatRasters")
		}
		x@ptr <- x@ptr$geometry(nlyrs, props, keeptime, keepunits, tags)
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
	function(x, crs="", extent=NULL) {
		dims <- dim(x)
		if (length(dims) < 3) {
			error("rast,array", "cannot handle an array with less than 3 dimensions")
		}
		if (length(dims) > 3) {
			if (length(dims) == 4) {
				if (dims[4] == 1) {
					x <- x[,,,1]
				} else {
					error("rast,array", "rast cannot handle an array with 4 dimensions (try 'sds')")
				}
			} else {
				error("rast,array", "cannot handle an array with more than 3 dimensions")
			}
		}
		r <- methods::new("SpatRaster")
		if (!is.null(extent)) {
			e <- as.vector(extent)
		} else {
			e <- c(0, dims[2], 0, dims[1])
		}
		crs <- character_crs(crs, "rast")
		r@ptr <- SpatRaster$new(dims, e, crs)
		values(r) <- x
		messages(r, "rast")
	}
)



setMethod("rast", signature(x="ANY"),
	function(x, ...) {
		if (inherits(x, "sf")) {
			out <- rast(ext(x), ...)
			if (is.null(list(...)$crs)) {
				sfi <- attr(x, "sf_column")
				crs(out, warn=FALSE) <- attr(x[[sfi]], "crs")$wkt
			}
		} else {
			out <- methods::as(x, "SpatRaster")
		}
		#g <- gc()
		out
	}
)


.rastFromXYZ <- function(xyz, digits=6, crs="", extent=NULL) {


	ln <- colnames(xyz)
	## xyz might not have colnames, or might have "" names
	if (is.null(ln)) ln <- rep("", ncol(xyz))
	if (any(nchar(ln) < 1)) ln <- make.names(ln)
	if (inherits(xyz, "data.frame")) {
		xyz <- as.matrix(xyz)
		xyz <- matrix(as.numeric(xyz), ncol=ncol(xyz), nrow=nrow(xyz))
	}
	x <- sort(unique(xyz[,1]))
	if (length(x) == 1) {
		error("rast", "cannot create a raster geometry from a single x coordinate")
	}
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
	if (length(y) == 1) {
		error("rast", "cannot create a raster geometry from a single y coordinate")
	}
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
		v <- try( matrix(NA, nrow=ncell(r), ncol= nlyr(r)) )
		if (inherits(v, "try-error")) {
			error(paste("cannot make matrix with ", ncell(r), " rows"))
		}
		v[cells, ] <- xyz[, -c(1:2)]
		values(r) <- v
	}

	if (!is.null(extent)) {
		r <- extend(r, extent)
		r <- crop(r, extent)
	}

	return(r)
}


.rastFromXYLZ <- function(x, digits=6, crs="", extent=NULL) {
	if (ncol(x) != 4) {
		error("rast", "a 'xylz' structure must have 4 columns")
	}
	nms <- names(x)
	names(x)[1:3] <- c("x", "y", "l")
	w <- stats::reshape(x, timevar="l", idvar=c("x", "y"), direction="wide")
	w <- rast(w, type="xyz", digits=digits, crs=crs, extent=extent)
	names(w) <- gsub(paste0(nms[4], "."), "", names(w))
	if (inherits(x[,3], "Date") || inherits(x[,3], "POSIXlt") || inherits(x[,3], "POSIXct")) {
		time(w) <- unique(x[,3])
		names(w) <- paste0(nms[4], ".", 1:nlyr(w))
	} else if (inherits(x[,3], "numeric")) {
		u <- unique(x[,3])
		if (all(u == trunc(u))) {
			su <- sort(u)
			if ((su[1] == 1) && (su[length(su)] == length(su))) {
				if (!all(su == u)) {
					w <- w[[order(u)]]
				}
			}
		}
	}
	w
}

	

setMethod("rast", signature(x="matrix"),
	function(x, type="", crs="", digits=6, extent=NULL) {
		stopifnot(prod(dim(x)) > 0)
		if (type == "xyz") {
			r <- .rastFromXYZ(x, crs=crs, digits=digits, extent=extent)
		} else if (type == "xylz") {
			r <- .rastFromXYLZ(x, crs=crs, digits=digits, extent=extent)
		} else if (type != "") {
			error("rast", 'argument type should be one of "", "xyz", or "xylz"')
		} else {
			if (is.null(extent)) {
				r <- rast(nrows=nrow(x), ncols=ncol(x), extent=ext(c(0, ncol(x), 0, nrow(x))), crs=crs)
			} else {
				r <- rast(nrows=nrow(x), ncols=ncol(x), crs=crs, extent=extent)
			}
			values(r) <- as.vector(t(x))
		}
		messages(r, "rast")
	}
)


setMethod("rast", signature(x="data.frame"),
	function(x, type="xyz", crs="", digits=6, extent=NULL) {
		if (type == "xyz") {
			.rastFromXYZ(x, crs=crs, digits=digits, extent=extent)
		} else if (type == "xylz") {
			r <- .rastFromXYLZ(x, crs=crs, digits=digits, extent=extent)
		} else if (type != "") {
			error("rast", 'argument type should be one of "", "xyz", or "xylz"')
		} else {
			rast(as.matrix(x), type=type, crs=crs, digits=digits, extent=extent)
		}
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


