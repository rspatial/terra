# Author: Robert J. Hijmans
# Date :  October 2017
# Version 1.0
# License GPL v3


setMethod("rast", signature(x="missing"),
	function(x, nrows=180, ncols=360, nlyrs=1, xmin=-180, xmax=180, ymin=-90, ymax=90, crs, extent, resolution, vals, ...) {

		if (missing(extent)) {	
			e <- c(xmin, xmax, ymin, ymax) 
		} else {
			e <- as.vector(extent)
		}
		if ((e[1] >= e[2]) || e[3] >= e[4]) {stop("invalid extent")}
		
		if (missing(crs)) {
			if (e[1] > -360.01 & e[2] < 360.01 & e[3] > -90.01 & e[4] < 90.01) {
				crs <- "+proj=longlat +datum=WGS84"
				#crs <- 'GEOGCS["WGS 84", DATUM["WGS_1984", SPHEROID["WGS 84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]]'
			} else {
				crs <- ""
			}
		} else {
			crs <- as.character(crs)
		}

		if (!missing(resolution)) {
			nrows <- max(1, round((e[4] - e[3]) / resolution))
			ncols <- max(1, round((e[2] - e[1]) / resolution))
		}

		r <- methods::new("SpatRaster")
		r@ptr <- SpatRaster$new(c(nrows, ncols, nlyrs), e, crs)
		r <- show_messages(r, "rast")
		if (!missing(vals)) {
			if (length(vals) == ncell(r)) {
				values(r) <- vals
			} else {
				values(r) <- rep_len(vals, ncell(r))
			}
		}
		r
	}
)

setMethod("rast", signature(x="list"),
	function(x, ...) {
		i <- sapply(x, function(i) inherits(i, "SpatRaster"))
		if (!any(i)) {
			stop("None of the elements of x are a SpatRaster")
		}
		if (!all(i)) {
			warning(paste(sum(!i), "out of", length(x), "elements of x are a SpatRaster"))
		}
		names(x) <- NULL
		do.call(c, x[i])
	}
)


setMethod("rast", signature(x="SpatExtent"),
	function(x, ...) {
		e <- as.vector(x)
		rast(xmin=e[1], xmax=e[2], ymin=e[3], ymax=e[4], ...)
	}
)

setMethod("rast", signature(x="SpatVector"),
	function(x, ...) {
		e <- ext(x)
		dx <- diff(as.vector(e)[1:2])
		dy <- diff(as.vector(e)[3:4])
		res <- min(dx, dy) / 100
		rast(ext(x), res=res, ...)
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
	function(x, subds=0, ...) {

		x <- trimws(x)
		x <- x[x!=""]
		if (length(x) == 0) {
			stop("provide a valid filename")
		}
		r <- methods::new("SpatRaster")
		f <- .fullFilename(x)
		#subds <- subds[1]
		if (is.character(subds)) { 
			r@ptr <- SpatRaster$new(f, -1, subds, "")		
		} else {
			r@ptr <- terra:::SpatRaster$new(f, subds-1, "", "")
		}
		if (r@ptr$getMessage() == "ncdf extent") {
			test <- try(r <- .ncdf_extent(r), silent=TRUE)
			if (inherits(test, "try-error")) {
				warning("GDAL did not find an extent. Cells not equally spaced?") 
			}
		}
		r <- show_messages(r, "rast")
		
		if (crs(r) == "") {
			if (isLonLat(r, perhaps=TRUE, warn=FALSE)) {
				crs(r) <- "+proj=longlat +datum=WGS84"
			}
		}
		r
	}
)


setMethod("rast", signature(x="SpatRaster"),
	function(x, nlyrs=nlyr(x), ...) {
		r <- methods::new("SpatRaster")
		r@ptr <- x@ptr$geometry(nlyrs, FALSE)
		if (length(list(...)) > 0) {
			warning("additional arguments are ignored")
		}
		show_messages(r, "rast")
	}
)


setMethod("rast", signature(x="SpatRasterDataset"),
	function(x, ...) {
		rast(x[1], ...)
	}
)



setMethod("rast", signature(x="array"),
	function(x, ...) {
		dims <- dim(x)
		if (length(dims) > 3) {
			stop("cannot handle an array with more than 3 dimensions")
		}
		r <- methods::new("SpatRaster")
		r@ptr <- SpatRaster$new(dims, c(0, dims[2], 0, dims[1]), "")
		values(r) <- x
		show_messages(r, "rast")
	}
)



setMethod("rast", signature(x="Raster"),
	function(x, ...) {
		methods::as(x, "SpatRaster")
	}
)


.rastFromXYZ <- function(xyz, digits=6, crs="", ...) {

	if (length(list(...))>0) warning("additional arguments ignored when x is a SpatRasterDataset")

	ln <- colnames(xyz)
	## xyz might not have colnames, or might have "" names
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
		stop("x cell sizes are not regular")
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
		stop("y cell sizes are not regular")
	}

	minx <- min(x) - 0.5 * rx
	maxx <- max(x) + 0.5 * rx
	miny <- min(y) - 0.5 * ry
	maxy <- max(y) + 0.5 * ry

	d <- dim(xyz)
	r <- rast(xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, crs=crs, nl=d[2]-2)
	res(r) <- c(rx, ry)
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
	function(x, type="", crs="", ...) {
		if (type == "xyz") {
			r <- .rastFromXYZ(x, ...)
		} else {
			r <- rast(nrows=nrow(x), ncols=ncol(x), extent=ext(c(0, ncol(x), 0, nrow(x))), crs=crs, ...)
			values(r) <- t(x)
		}
		show_messages(r, "rast")
	}
)


setMethod("NAflag<-", signature(x="SpatRaster"), 
	function(x, ..., value)  {
		value <- as.numeric(value)
		if (!(x@ptr$setNAflag(value))) {
			stop("cannot set these values")
		}
		x
	}
)

setMethod("NAflag", signature(x="SpatRaster"), 
	function(x, ...)  {
		x@ptr$getNAflag()
	}
)

