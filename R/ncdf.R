
.ncdf_extent <- function(x) {

#	stopifnot(requireNamespace("ncdf4"))
	s <- sources(x)[1,1]

	ss <- unlist(strsplit(s, "\""))
	zvar <- gsub(":", "", ss[3])
	dims <- 1:3

	nc <- ncdf4::nc_open(ss[2], readunlim=FALSE, suppress_dimvals = TRUE)
	on.exit( ncdf4::nc_close(nc) )		

	ncols <- nc$var[[zvar]]$dim[[dims[1]]]$len
	nrows <- nc$var[[zvar]]$dim[[dims[2]]]$len
	if (!(ncol(x) == ncols) & (nrow(x) == nrows)) {
		warning("GDAL did not find an extent. Cells not equally spaced?") 
		return(x)
	}
	xx <- try(ncdf4::ncvar_get(nc, nc$var[[zvar]]$dim[[dims[1]]]$name), silent = TRUE)
	if (inherits(xx, "try-error")) {
	  xx <- seq_len(nc$var[[zvar]]$dim[[dims[1]]]$len)
	}
	rs <- xx[-length(xx)] - xx[-1]
	if (! isTRUE ( all.equal( min(rs), max(rs), tolerance=0.025, scale= abs(min(rs))) ) ) {
		warning("cells are not equally spaced; extent is not defined") 
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
		warning("cells are not equally spaced; extent is not defined") 
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
	function(x, filename, varname, longname="", units="", overwrite=FALSE, datatype="double", NAflag=-9999, ...) {
		units(x) <- units
		varnames(x) <- varname
		longnames(x) <- longname
		x <- sds(x)
		writeCDF(x, filename=filename, overwrite=overwrite, datatype=datatype, NAflag=NAflag, ...)
	}
)


setMethod("writeCDF", signature(x="SpatDataSet"), 
	function(x, filename, overwrite=FALSE, datatype="double", NAflag=-9999, ...) {

		force_v4 <- list(...)$force_v4
		if (is.null(force_v4)) force_v4 <- FALSE

		filename <- trimws(filename)
		stopifnot(filename != "")
		if (file.exists(filename) & !overwrite) {
			stop("file exists, use overwrite=TRUE to overwrite it")
		}

		n <- length(x)
#		longnames <- longnames(x)
		vars <- varnames(x)
		vars[vars == ""] <- (1:n)[vars == ""] 
#		units <- units(x)
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
				unit <- units(y)[1]
				lvar <- longnames(y)[1]
				if (y@ptr$hasTime) {
					zv <- time(y)
					zatt <- list("units=seconds since 1970-1-1 00:00:00")		
					zunit <- "seconds"
					zname <- "time"
				} else {
					zv <- 1:nlyr(y)
					zatt <- list("units=unknown")		
					zunit <- "unknown"
					zname <- "layer"
				} 
				zdim <- ncdf4::ncdim_def(zname, zunit, zv, unlim=TRUE )

				ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units, list(xdim, ydim, zdim), NAflag, lvar, prec = dtype[i], ...)
			} else {			
				ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units, list(xdim, ydim), NAflag, lvar, prec = dtype[i], ...)
			}
		}
		
		ncobj <- ncdf4::nc_create(filename, ncvars, force_v4=force_v4)
		#on.exit( ncdf4::nc_close(ncobj) )
		
		# writing all at once. need to chunk 
		nc <- ncol(x)
		nr <- nrow(x)
		for (i in 1:n) {
			y = x[i]
			if (nl[i] > 1) {
				v <- values(y)
				v <- array(v, c(nr, nc, nl[i]))

				ncdf4::ncvar_put(ncobj, ncvars[[i]], values(y), start=c(1,1,1), count=dim(y)[c(2,1,3)])
			} else {
				v <- t(as.matrix(y, TRUE))
				ncdf4::ncvar_put(ncobj, ncvars[[i]], v, start=c(1,1), count=dim(y)[2:1])			
			}
		}

		#crsname <- "latitude_longitude"
		#ncvars[[7]] <- ncdf4::ncvar_def(crsname, "1", NULL, NULL, crsname, prec = "char")
		#nc <- nc_create(filename, ncvars)
		#prj <- crs(x[1])
		#if (!is.na(prj)) {
		#	ncdf4::ncatt_put(ncobj, "crs", "wkt", prj, prec="text")
		#	ncdf4::ncatt_put(ncobj, varname, "grid_mapping", "crs")
		#	ncdf4::ncatt_put(ncobj, varname, "wkt", prj, prec="text")
		#}
		
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

