
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


.getNetCDFDType <- function(dtype) {
	if (!(dtype %in% c('LOG1S', 'INT1S', 'INT2S', 'INT4S', 'INT1U', 'INT2U', 'FLT4S', 'FLT8S'))) {
		stop('not a valid data type')
	}
	type <- substr(dtype, 1, 3)
	size <- as.integer(substr(dtype, 4, 4)) * 8
	signed <- substr(dtype, 5, 5) == "S"
	
	if (size == 8) {
		if (!signed) {
			return("char") #8-bit characters intended for representing text.
		} else {
			return("byte")
		}
	} else if (type == 'INT') {
		if (!signed) {
			warning('netcdf only stores signed integers')
		}
		if (size == 16) { 
			return( "short" ) 
		} else if (size == 32 ) { 
			return( "integer" ) 
		} else {
			return ( "double" )		
		}
	} else {
		if (size == 32) { 
			return( "float" ) 
		} else {  
			return ( "double" )  
		}
	}
}




.writeNCDF <- function(x, filename, datatype='FLT4S', overwrite=FALSE, att, varname, varunit, varatt, longname, xname, yname, zname, zunit, zatt, NAflag=-9999, compr=6, force_v4=TRUE, ...) {
		
	filename <- trim(filename)
	stopifnot(filename != "")
	if (file.exists(filename) & !overwrite) {
		stop("file exists, use overwrite=TRUE to overwrite it")
	}
	ncdatatype <- .getNetCDFDType(datatype)
	nl <- nlyr(x)
	
	if (isLonLat(x, perhaps=TRUE, warn=FALSE)) {
		if (missing(xname)) xname = 'longitude'
		if (missing(yname)) yname = 'latitude'
		xunit = 'degrees_east'
		yunit = 'degrees_north'
	} else {
		if (missing(xname)) xname = 'easting'
		if (missing(yname)) yname = 'northing'
		xunit = 'meter' # probably
		yunit = 'meter' # probably
	}
	
	if (missing(varname))  {
		varname <- gsub("_1$", "", names(x)[1])
		varname <- gsub(".1$", "", names(x)[1])
	}
	if (missing(varunit))  varunit <- ""
	if (missing(longname))  longname <- ""

	ht <- x@ptr$hasTime
	zv <- x@ptr$time
	zatt <- list('units=seconds since 1970-1-1 00:00:00')		
	zunit <- 'seconds'
	zname <- "time"
	
	xdim <- ncdf4::ncdim_def( xname, xunit, xFromCol(x, 1:ncol(x)) )
	ydim <- ncdf4::ncdim_def( yname, yunit, yFromRow(x, 1:nrow(x)) )
	zdim <- ncdf4::ncdim_def( zname, zunit, zv, unlim=TRUE )
	vardef <- ncdf4::ncvar_def( varname, varunit, list(xdim, ydim, zdim), NAflag, longname, prec = ncdatatype, ... )
	crsdef <- ncdf4::ncvar_def("crs", "", list(), NULL, prec="integer")
	defs <- list(crsdef, vardef)

	nc <- ncdf4::nc_create(filename, defs, force_v4=force_v4)
	on.exit( ncdf4::nc_close(nc) )		

	prj <- crs(x)
	if (!is.na(prj)) {
		ncdf4::ncatt_put(nc, "crs", "wkt", prj, prec='text')
		ncdf4::ncatt_put(nc, varname, "grid_mapping", "crs")
		ncdf4::ncatt_put(nc, varname, "wkt", prj, prec='text')
	}
	ncdf4::ncatt_put(nc, 0, 'Conventions', 'CF-1.4', prec='text')

	pkgversion <- drop(read.dcf(file=system.file("DESCRIPTION", package='terra'), fields=c("Version")))
	ncdf4::ncatt_put(nc, 0, 'created_by', paste('R, packages ncdf4 and terra (version ', pkgversion, ')', sep=''), prec='text')
	ncdf4::ncatt_put(nc, 0, 'date', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), prec='text')

#start loop
	lstart <- 1
	start <- 1
	lend <- nl
	v <- values(x)
	nrows <- nrow(x)
	ncols <- ncol(x)
	v <- array(v, c(nrows, ncols, nl))
	try ( ncdf4::ncvar_put(nc, varname, v, start=c(1, start, lstart), count=c(ncols, nrows, lend) ) )
#end loop	

	rast(filename)
}

