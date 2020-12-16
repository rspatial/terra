
.ncdf_extent <- function(x) {

	if (!("ncdf4" %in% rownames(installed.packages()))) {
		warn("rast", "GDAL did not find an extent. Cells not equally spaced?") 
		warn("rast", "installing the ncdf4 package may help")
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


setMethod("writeCDF", signature(x="SpatRaster"), 
	function(x, filename, varname, longname="", unit="", overwrite=FALSE, datatype="double", NAflag=-9999, ...) {
		filename <- trimws(filename)
		stopifnot(filename != "")
		if (missing(varname)) {
			varname <- tools::file_path_sans_ext(basename(filename))
		}
		varnames(x) <- varname
		longnames(x) <- longname
		units(x) <- unit
		x <- sds(x)
		writeCDF(x, filename=filename, overwrite=overwrite, datatype=datatype, NAflag=NAflag, ...)
	}
)


setMethod("writeCDF", signature(x="SpatRasterDataset"), 
	function(x, filename, overwrite=FALSE, datatype="double", NAflag=-9999, zname="time", ...) {
		filename <- trimws(filename)
		stopifnot(filename != "")

		dots <- list(...)
		force_v4 <- if (is.null(dots$force_v4)) { TRUE } else {dots$force_v4}
		verbose <- if (is.null(dots$verbose)) { FALSE } else  { dots$verbose }
		filename <- trimws(filename)
		stopifnot(filename != "")
		if (file.exists(filename) & !overwrite) {
			error("writeCDF", "file exists, use 'overwrite=TRUE' to overwrite it")
		}

		n <- length(x)
		lvar <- longnames(x)
		vars <- varnames(x)
		vars[vars == ""] <- paste0("var_", (1:n)[vars == ""])
		units <- units(x)
		nl <- nlyr(x)
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

		dtype <- rep_len(datatype, n)
		ncvars <- list()
		for (i in 1:n) {
			y <- x[i]
			if (nl[i] > 1) {	
				lvar <- longnames(y)[1]
				if (y@ptr$hasTime) {
					zv <- y@ptr$time
					zunit <- "seconds since 1970-1-1 00:00:00"
					zname <- "time"
					cal <- "standard"
				} else {
					zv <- 1:nlyr(y)
					zunit <- "unknown"
					cal <- NA
				} 
				zdim <- ncdf4::ncdim_def(zname, zunit, zv, unlim=TRUE, calendar=cal)
				ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units[i], list(xdim, ydim, zdim), NAflag, lvar, prec = dtype[i], ...)
			} else {			
				ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units[i], list(xdim, ydim), NAflag, lvar, prec = dtype[i], ...)
			}
		}
		ncvars[[n+1]] <- ncdf4::ncvar_def("crs", "", list(), NULL, prec="integer")
		
		ncobj <- ncdf4::nc_create(filename, ncvars, force_v4=force_v4, verbose=verbose)

		prj <- crs(x[1])
		prj <- gsub("\n", "", prj)
		if (prj != "") {
			ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "spatial_ref", prj, prec="text")
			ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "proj4", terra:::.proj4(x[1]), prec='text')
		}
		e <- ext(x)
		rs <- res(x)
		gt <- paste(trimws(formatC(as.vector(c(e$xmin, rs[1], 0, e$ymax, 0, -1 * rs[2])), 22)), collapse=" ")
		ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "GeoTransform", gt, prec="text")

		nc <- ncol(x)
		nr <- nrow(x)
		nl <- nlyr(x)

		opt <- spatOptions("", TRUE, list())

		for (i in 1:n) {
			y = x[i]
			b <- y@ptr$getBlockSize(4, opt$memfrac)
			if (nl[i] > 1) {
				for (j in 1:b$n) {
					d <- readValues(y, b$row[j]+1, b$nrows[j], 1, nc, FALSE, FALSE)
					d <- array(d, c(nc, nr, nl))		
					ncdf4::ncvar_put(ncobj, ncvars[[i]], d, start=c(1, b$row[j]+1, 1), count=c(nc, b$nrows[j], nl))
				}
			} else {
				for (j in 1:b$n) {
					d <- readValues(y, b$row[j]+1, b$nrows[j], 1, nc, FALSE, FALSE)
					d <- matrix(d, ncol=b$nrows[j])
					ncdf4::ncvar_put(ncobj, ncvars[[i]], d, start=c(1, b$row[j]+1), count=c(nc, b$nrows[j]))
				}
			}
			if (prj != "") {
				ncdf4::ncatt_put(ncobj, ncvars[[i]], "grid_mapping", "crs", prec="text")
			}
		}

		
		ncdf4::ncatt_put(ncobj, 0, "Conventions", "CF-1.4", prec="text")
		pkgversion <- drop(read.dcf(file=system.file("DESCRIPTION", package="terra"), fields=c("Version")))
		ncdf4::ncatt_put(ncobj, 0, "created_by", paste("R, packages ncdf4 and terra (version ", pkgversion, ")", sep=""), prec="text")
		ncdf4::ncatt_put(ncobj, 0, "date", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), prec="text")
		ncdf4::nc_close(ncobj)

		invisible(rast(filename))
	}
)



.vectCDF <- function(filename, varname, polygons=FALSE) {
# read (irregular) raster netcdf as points or polygons
# not to be confused with vector netcdf format
	

}

