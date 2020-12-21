
.ncdf_extent <- function(x) {

	if (!("ncdf4" %in% rownames(utils::installed.packages()))) {
		warn("rast", "GDAL did not find an extent. installing the ncdf4 package may help")
		return(x)
	}
	s <- sources(x)[1,1]

	ss <- unlist(strsplit(s, "\""))
	zvar <- gsub(":", "", ss[3])
	dims <- 1:3

	nc <- ncdf4::nc_open(ss[2], readunlim=FALSE, suppress_dimvals = TRUE)
	on.exit( ncdf4::nc_close(nc) )		

	ncols <- nc$var[[zvar]]$dim[[dims[1]]]$len
	nrows <- nc$var[[zvar]]$dim[[dims[2]]]$len
	if (!(ncol(x) == ncols) & (nrow(x) == nrows)) {
		warn("rast", "GDAL did not find an extent. Cells not equally spaced?") 
		return(x)
	}
	xx <- try(ncdf4::ncvar_get(nc, nc$var[[zvar]]$dim[[dims[1]]]$name), silent = TRUE)
	if (inherits(xx, "try-error")) {
	  xx <- seq_len(nc$var[[zvar]]$dim[[dims[1]]]$len)
	}
	rs <- xx[-length(xx)] - xx[-1]
	if (! isTRUE ( all.equal( min(rs), max(rs), tolerance=0.025, scale= abs(min(rs))) ) ) {
		warn("rast", "cells are not equally spaced; extent is not defined") 
		return(x)
	}
	xrange <- c(min(xx), max(xx))
	resx <- (xrange[2] - xrange[1]) / (ncols-1)

	yy <- try(ncdf4::ncvar_get(nc, nc$var[[zvar]]$dim[[dims[2]]]$name), silent = TRUE)
	if (inherits(yy, "try-error")) {
	  yy <- seq_len(nc$var[[zvar]]$dim[[dims[2]]]$len)
	}
	
	rs <- yy[-length(yy)] - yy[-1]
	if (! isTRUE ( all.equal( min(rs), max(rs), tolerance=0.025, scale= abs(min(rs))) ) ) {
		warn("rast", "cells are not equally spaced; extent is not defined") 
		return(x)
	}
	yrange <- c(min(yy), max(yy))
	resy <- (yrange[2] - yrange[1]) / (nrows-1)

	xrange[1] <- xrange[1] - 0.5 * resx
	xrange[2] <- xrange[2] + 0.5 * resx
	yrange[1] <- yrange[1] - 0.5 * resy
	yrange[2] <- yrange[2] + 0.5 * resy
	ext(x) <- ext(xrange[1], xrange[2], yrange[1], yrange[2])
	return(x)
}



.write_cdf <- function(x, filename, overwrite=FALSE, zname="time", missval=-9999, prec="float", compression=NA, ...) {

	dots <- list(...)
	force_v4 <- if (is.null(dots$force_v4)) { TRUE } else {dots$force_v4}
	verbose  <- if (is.null(dots$verbose)) { FALSE } else {dots$verbose}

	n <- length(x)
	y <- x[1]
	if (isLonLat(y, perhaps=TRUE, warn=FALSE)) {
		xname = "longitude"
		yname = "latitude"
		xunit = "degrees_east"
		yunit = "degrees_north"
	} else {
		xname = "easting"
		yname = "northing"
		xunit = "meter" # probably
		yunit = "meter" # probably
	}
	xdim <- ncdf4::ncdim_def( xname, xunit, xFromCol(y, 1:ncol(y)) )
	ydim <- ncdf4::ncdim_def( yname, yunit, yFromRow(y, 1:nrow(y)) )
	vars <- varnames(x)
	vars[vars == ""] <- paste0("var_", (1:n)[vars == ""])
	vars <- make.unique(vars)

	lvar <- longnames(x)
	units <- units(x)
	zname <- rep_len(zname, n)

	prec <- rep_len(prec, n)
	compression <- rep_len(compression, n)
	nc <- ncol(x)
	nr <- nrow(x)
	nl <- nlyr(x)
	ncvars <- list()
	cal <- NA
	for (i in 1:n) {
		if (nl[i] > 1) {	
			y <- x[i]
			if (y@ptr$hasTime) {
				zv <- y@ptr$time
				tstep <- y@ptr$timestep 
				cal <- "standard"
				if (tstep == "seconds") {
					zunit <- "seconds since 1970-1-1 00:00:00"
					cal <- "standard"
				} else if (tstep == "days") {
					zunit <- "days since 1970-1-1"
					cal <- "standard"
				} else {
					zunit <- "unknown"					
				}
			} else {
				zv <- 1:nlyr(y)
				zunit <- "unknown"
			} 
			zdim <- ncdf4::ncdim_def(zname[i], zunit, zv, unlim=FALSE, create_dimvar=TRUE, calendar=cal)
			ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units[i], list(xdim, ydim, zdim), missval, lvar[i], prec = prec[i], compression=compression[i],...)
		} else {			
			ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units[i], list(xdim, ydim), missval, lvar[i], prec = prec[i], compression=compression[i], ...)
		}
	}

	ncvars[[n+1]] <- ncdf4::ncvar_def("crs", "", list(), NULL, prec="integer")
		
	ncobj <- ncdf4::nc_create(filename, ncvars, force_v4=force_v4, verbose=verbose)
	on.exit(ncdf4::nc_close(ncobj))

	prj <- crs(x[1])
	prj <- gsub("\n", "", prj)
	if (prj != "") {
		ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "spatial_ref", prj, prec="text")
		ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "proj4", .proj4(x[1]), prec='text')
	}
	e <- ext(x)
	rs <- res(x)
	gt <- paste(trimws(formatC(as.vector(c(e$xmin, rs[1], 0, e$ymax, 0, -1 * rs[2])), 22)), collapse=" ")
	ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "GeoTransform", gt, prec="text")


	opt <- spatOptions("", TRUE, list())

	for (i in 1:n) {
		y = x[i]
		readStart(y)
		b <- y@ptr$getBlockSize(4, opt$memfrac)
		if (nl[i] > 1) {
			for (j in 1:b$n) {
				d <- readValues(y, b$row[j]+1, b$nrows[j], 1, nc, FALSE, FALSE)
				d <- array(d, c(nc, b$nrows[j], nl[i]))		
				ncdf4::ncvar_put(ncobj, ncvars[[i]], d, start=c(1, b$row[j]+1, 1), count=c(nc, b$nrows[j], nl[i]))
			}
		} else {
			for (j in 1:b$n) {
				d <- readValues(y, b$row[j]+1, b$nrows[j], 1, nc, FALSE, FALSE)
				d <- matrix(d, ncol=b$nrows[j])
				ncdf4::ncvar_put(ncobj, ncvars[[i]], d, start=c(1, b$row[j]+1), count=c(nc, b$nrows[j]))
			}
		}
		readStop(y)
		if (prj != "") {
			ncdf4::ncatt_put(ncobj, ncvars[[i]], "grid_mapping", "crs", prec="text")
		}
	}
		
	ncdf4::ncatt_put(ncobj, 0, "Conventions", "CF-1.4", prec="text")
	pkgversion <- drop(read.dcf(file=system.file("DESCRIPTION", package="terra"), fields=c("Version")))
	ncdf4::ncatt_put(ncobj, 0, "created_by", paste("R, packages ncdf4 and terra (version ", pkgversion, ")", sep=""), prec="text")
	ncdf4::ncatt_put(ncobj, 0, "date", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), prec="text")
	
	TRUE
}



setMethod("writeCDF", signature(x="SpatRaster"), 
	function(x, filename, varname, longname="", unit="", ...) {
		filename <- trimws(filename)
		stopifnot(filename != "")
		if (missing(varname)) {
			varname <- tools::file_path_sans_ext(basename(filename))
		}
		varnames(x) <- varname
		longnames(x) <- longname
		units(x) <- unit
		x <- sds(x)
		invisible( writeCDF(x, filename=filename, ...) )
	}
)


setMethod("writeCDF", signature(x="SpatRasterDataset"), 
	function(x, filename, overwrite=FALSE, zname="time", missval=-9999, prec="float", compression=NA, ...) {
		filename <- trimws(filename)
		stopifnot(filename != "")
		xt  <- file_ext(filename)
		if (!(xt %in% c(".nc", ".cdf"))) {
			warn("writeCDF", "for better results use file extension '.nc' or '.cdf'\nsee: https://stackoverflow.com/a/65398262/635245")
		}
		if (file.exists(filename) & !overwrite) {
			error("writeCDF", "file exists, use 'overwrite=TRUE' to overwrite it")
		}
		ok <- .write_cdf(x, filename, zname=zname, missval=missval, prec=prec, compression=compression, ...)
		if (ok) {
			if (length(x) > 1) {
				out <- sds(filename)
			} else {
				out <- rast(filename)
			}
			invisible(out)
		} else {
			error("writeCDF", "?")
		}
	}
)




#.vectCDF <- function(filename, varname, polygons=FALSE) {
  # read (irregular) raster netcdf as points or polygons
  # not to be confused with vector netcdf format
#}

