
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
	function(x, filename, overwrite=FALSE, datatype="double", NAflag=-9999, ...) {
		x <- sds(x)
		writeCDF(x, filename=filename, overwrite=overwrite, datatype=datatype, NAflag=NAflag, ...)
	}
)


setMethod("writeCDF", signature(x="SpatDataSet"), 
	function(x, filename, overwrite=FALSE, datatype="double", NAflag=-9999, ...) {

		force_v4=TRUE

		filename <- trimws(filename)
		stopifnot(filename != "")
		if (file.exists(filename) & !overwrite) {
			stop("file exists, use overwrite=TRUE to overwrite it")
		}

		# loop over subdatasets 
		# for now:
		x <- x[1] 
		longname <- longnames(x)[1]
		varname <- varnames(x)[1]
		if (varname == "") varname <- "data"
		unit <- units(x)[1]
		
		
		if (isLonLat(x, perhaps=TRUE, warn=FALSE)) {
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
		xdim <- ncdf4::ncdim_def( xname, xunit, xFromCol(x, 1:ncol(x)) )
		ydim <- ncdf4::ncdim_def( yname, yunit, yFromRow(x, 1:nrow(x)) )

		nl <- nlyr(x)
	
		if (nl > 1) {
			if (x@ptr$hasTime) {
				zv <- x@ptr$time
				zatt <- list("units=seconds since 1970-1-1 00:00:00")		
				zunit <- "seconds"
				zname <- "time"
			} else {
				zv <- 1:nlyr(x)
				zatt <- list("units=unknown")		
				zunit <- "unknown"
				zname <- "layer"
			}
			zdim <- ncdf4::ncdim_def( zname, zunit, zv, unlim=TRUE )
			vardef <- ncdf4::ncvar_def( varname, unit, list(xdim, ydim, zdim), NAflag, longname, prec = datatype, ... )
		} else {
			vardef <- ncdf4::ncvar_def( varname, unit, list(xdim, ydim), NAflag, longname, prec = datatype, ... )		
		}
		
		crsdef <- ncdf4::ncvar_def("crs", "", list(), NULL, prec="integer")
		defs <- list(crsdef, vardef)

		nc <- ncdf4::nc_create(filename, defs, force_v4=force_v4)
		on.exit( ncdf4::nc_close(nc) )		

		prj <- crs(x)
		if (!is.na(prj)) {
			ncdf4::ncatt_put(nc, "crs", "wkt", prj, prec="text")
			ncdf4::ncatt_put(nc, varname, "grid_mapping", "crs")
			ncdf4::ncatt_put(nc, varname, "wkt", prj, prec="text")
		}
		
		ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.4", prec="text")
		pkgversion <- drop(read.dcf(file=system.file("DESCRIPTION", package="terra"), fields=c("Version")))
		ncdf4::ncatt_put(nc, 0, "created_by", paste("R, packages ncdf4 and terra (version ", pkgversion, ")", sep=""), prec="text")
		ncdf4::ncatt_put(nc, 0, "date", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), prec="text")

	#start loop for writing
	# for now:
		nrows <- nrow(x)
		ncols <- ncol(x)
		start <- 1
		v <- values(x)

		if (nl > 1) {
			lstart <- 1
			lend <- nl
			v <- array(v, c(nrows, ncols, nl))
			try ( ncdf4::ncvar_put(nc, varname, v, start=c(1, start, lstart), count=c(ncols, nrows, lend) ) )
		} else {
			v <- array(v, c(nrows, ncols, nl))
			try ( ncdf4::ncvar_put(nc, varname, v, start=c(1, start), count=c(ncols, nrows) ) )
		}

	#end loop	

		invisible(rast(filename))
	}
)



.vectCDF <- function(filename, varname, polygons=FALSE) {
# read (irregular) raster netcdf as points or polygons
# not to be confused with vector netcdf format
	

}

