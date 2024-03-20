
.ncdf_extent <- function(x, f) {

	if (!("ncdf4" %in% rownames(utils::installed.packages()))) {
		warn("rast", "GDAL did not find an extent. installing the ncdf4 package may help")
		return(x)
	}
	fname <- f
	zvar <- varnames(x)[1]

	dims <- 1:3

	nc <- ncdf4::nc_open(fname, readunlim=FALSE, suppress_dimvals = TRUE)
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



.write_cdf <- function(x, filename, overwrite=FALSE, zname="time", atts="", gridmap="", prec="float", compression=NA, missval, force_v4=TRUE, verbose=FALSE, ...) {

	n <- length(x)
	y <- x[1]
	if (is.lonlat(y, perhaps=TRUE, warn=FALSE)) {
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

	valid_prec <- c("short", "integer", "float", "double", "byte")
	if (!all(prec %in% valid_prec)) {
		error("writeCDF", paste("prec must be one of:", paste(valid_prec, collapse=", ")))
	}
	prec <- rep_len(prec, n)
	if (missing(missval)) {
		miss_vals <- c(-32768, -2147483647, -1.175494e38, -1.7976931348623157e308, 255)
		missval <- miss_vals[match(prec, valid_prec)]
	} else {
		missval <- rep_len(missval, n)
	}
	compression <- compression[1]
	nc <- ncol(x)
	nr <- nrow(x)
	nl <- nlyr(x)
	ncvars <- list()
	cal <- NA
	for (i in 1:n) {
		if ((nl[i] > 1) || (x[i]@ptr$hasTime)) {
			y <- x[i]
			if (y@ptr$hasTime) {
				zv <- y@ptr$time
				tstep <- y@ptr$timestep
				cal <- "standard"
				if (tstep == "seconds") {
					zunit <- "seconds since 1970-1-1 00:00:00"
				} else if (tstep == "days") {
					zunit <- "days since 1970-1-1"
					zv <- zv / (24 * 3600)
				} else if (tstep == "months") {
					zunit <- "months"
					zv <- time(y)
				} else if (tstep == "yearmonths") {
					zunit <- "months since 1970"
					tm <- time(y) - 1970
					yr <- tm %/% 1
					zv <- (yr*12) + round(12 * (tm %% 1))
				} else if (tstep == "years") {
					zunit <- "years since 1970"
					zv <- time(y) - 1970
				} else {
					zunit <- "unknown"
				}
			} else {
				zv <- 1:nlyr(y)
				zunit <- "unknown"
			}
			zdim <- ncdf4::ncdim_def(zname[i], zunit, zv, unlim=FALSE, create_dimvar=TRUE, calendar=cal)
			ncvars[[i]] <- ncdf4::ncvar_def(vars[i], units[i], list(xdim, ydim, zdim), missval[i], lvar[i], prec = prec[i], compression=compression,...)
		} else {

			ncvars[[i]] <- ncdf4::ncvar_def(name=vars[i], units=units[i], dim=list(xdim, ydim), missval=missval[i], longname=lvar[i], prec = prec[i], compression=compression, ...)
		}
	}

	ncvars[[n+1]] <- ncdf4::ncvar_def("crs", "", list(), NULL, prec="integer")
	
	ncobj <- ncdf4::nc_create(filename, ncvars, force_v4=force_v4, verbose=verbose)
	on.exit(ncdf4::nc_close(ncobj))

	haveprj <- FALSE
	prj <- crs(x[1])
	prj <- gsub("\n", "", prj)
	if (prj != "") {
		haveprj <- TRUE
		ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "crs_wkt", prj, prec="text")
		# need for older gdal?
		ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "spatial_ref", prj, prec="text")
		prj <- .proj4(x[1])
		if (prj != "") {
			ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "proj4", prj, prec='text')
		}
		prj <- crs(x[1], describe=TRUE)[1,3]
		if (!is.na(prj)) {
			ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "epsg_code", prj, prec='text')
		}
	}
	gridmap <- grep("=", gridmap, value=TRUE)
	if (length(gridmap)>0) {
		gridmap <- strsplit(gridmap, "=")
		for (i in 1:length(gridmap)) {		
			ncdf4::ncatt_put(ncobj, ncvars[[n+1]], gridmap[[i]][1], gridmap[[i]][2], prec="text")
		}
		haveprj <- TRUE
	}


	e <- ext(x)
	rs <- res(x)
	gt <- paste(trimws(formatC(as.vector(c(e$xmin, rs[1], 0, e$ymax, 0, -1 * rs[2])), 22)), collapse=" ")
	ncdf4::ncatt_put(ncobj, ncvars[[n+1]], "geotransform", gt, prec="text")

	opt <- spatOptions()
	for (i in 1:n) {
		y = x[i]
		readStart(y)
		b <- blocks(y, 4)
		if (length(ncvars[[1]]$dim) == 3) {
			for (j in 1:b$n) {
				d <- readValues(y, b$row[j], b$nrows[j], 1, nc, FALSE, FALSE)
				d[is.nan(d)] <- NA
				d <- array(d, c(nc, b$nrows[j], nl[i]))
				ncdf4::ncvar_put(ncobj, ncvars[[i]], d, start=c(1, b$row[j], 1), count=c(nc, b$nrows[j], nl[i]))
			}
		} else {
			for (j in 1:b$n) {
				d <- readValues(y, b$row[j], b$nrows[j], 1, nc, FALSE, FALSE)
				d[is.nan(d)] <- NA
				d <- matrix(d, ncol=b$nrows[j])
				ncdf4::ncvar_put(ncobj, ncvars[[i]], d, start=c(1, b$row[j]), count=c(nc, b$nrows[j]))
			}
		}
		readStop(y)
		if (haveprj) {
			ncdf4::ncatt_put(ncobj, ncvars[[i]], "grid_mapping", "crs", prec="text")
		}
	}

	ncdf4::ncatt_put(ncobj, 0, "Conventions", "CF-1.4", prec="text")
	pkgversion <- drop(read.dcf(file=system.file("DESCRIPTION", package="terra"), fields=c("Version")))
	ncdf4::ncatt_put(ncobj, 0, "created_by", paste("R packages ncdf4 and terra (version ", pkgversion, ")", sep=""), prec="text")
	ncdf4::ncatt_put(ncobj, 0, "date", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), prec="text")

	atts <- grep("=", atts, value=TRUE)
	if (length(atts) > 0) {
		atts <- strsplit(atts, "=")
		for (i in 1:length(atts)) {
			ncdf4::ncatt_put(ncobj, 0, atts[[i]][1], atts[[i]][2], prec="text")
		}
	}
	TRUE
}



setMethod("writeCDF", signature(x="SpatRaster"),
	function(x, filename, varname, longname="", unit="", split=FALSE, ...) {
		filename <- trimws(filename)
		if (length(filename) > 1) {
			if (length(filename) != nlyr(x)) {
				stop("either provide a single filename, or the same number as nlyr(x)")
			}
		}
		stopifnot(filename != "")
		if (split) {
			y <- sds(as.list(x))
			names(y) <- names(x)
			if (!missing(varname)) {
				varnames(y) <- varname
			}
			if (!missing(longname)) {
				longnames(y) <- longname
			}
			if (!missing(unit)) {
				units(y) <- unit
			} else {
				units(y) <- units(x)
			}
			invisible( writeCDF(y, filename=filename, ...) )
		} else {
			if (missing(varname)) {
				varname <- tools::file_path_sans_ext(basename(filename))
			}
			varnames(x) <- varname
			longnames(x) <- longname
			units(x) <- unit
			x <- sds(x)
			invisible( writeCDF(x, filename=filename, ...) )
		}
	}
)


setMethod("writeCDF", signature(x="SpatRasterDataset"),
	function(x, filename, overwrite=FALSE, zname="time", atts="", gridmap="", prec="float", compression=NA, missval, ...) {
		filename <- trimws(filename)
		stopifnot(filename != "")
		xt  <- tools::file_ext(filename)
		if (any(!(xt %in% c("nc", "cdf")))) {
			warn("writeCDF", "for better results use file extension '.nc' or '.cdf'\nsee: https://stackoverflow.com/a/65398262/635245")
		}
		if (file.exists(filename) & !overwrite) {
			error("writeCDF", "file exists, use 'overwrite=TRUE' to overwrite it")
		}
		ok <- .write_cdf(x, filename, zname=zname, atts=atts, gridmap=gridmap, prec=prec, compression=compression, missval=missval, ...)
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



.varName <- function(nc, varname="", warn=TRUE) {
	n <- nc$nvars
	dims <- vars <- vector(length=n)
	if (n > 0) {
		for (i in 1:n) {
			vars[i] <- nc$var[[i]]$name
			dims[i] <- nc$var[[i]]$ndims
		}
		vars <- vars[dims > 1]
		dims <- dims[dims > 1]
	}

	if (varname=='') {
		nv <- length(vars)
		if (nv == 0) {
			return('z')
		}

		if (nv  == 1) {
			varname <- vars
		} else {
			varname <- vars[which.max(dims)]
			if (warn) {
				if (sum(dims == max(dims)) > 1) {
					vars <- vars[dims==max(dims)]
					warning('varname used is: ', varname, '\nIf that is not correct, you can set it to one of: ', paste(vars, collapse=", ") )
				}
			}
		}
	}

	zvar <- which(varname == vars)
	if (length(zvar) == 0) {
		stop('varname: ', varname, ' does not exist in the file. Select one from:\n', paste(vars, collapse=", ") )
	}
	return(varname)
}



.getCRSfromGridMap4 <- function(g) {
	if (!is.null(g$epsg_code)) {
		crs <- g$epsg_code
		if (!grep("EPSG:", crs, ignore.case=TRUE)) {
			crs <- paste0("epsg:", crs)
		}
		return(crs)
	}

	sp <- g$standard_parallel
	if (length(sp) > 1) {
		g$standard_parallel1 <- sp[1]
		g$standard_parallel2 <- sp[2]
		g$standard_parallel <- NULL
	}

	vals <- sapply(g, function(i) i[1])
	vars <- names(vals)
	if (any(vars %in% c("proj4", "crs_wkt", "spatial_ref"))) {
		crs=vals[vars %in% c("proj4", "crs_wkt", "spatial_ref")][1]
		return(crs)
	}
# based on info at
# http://trac.osgeo.org/gdal/wiki/NetCDF_ProjectionTestingStatus
# accessed 7 October 2012
	prj <- matrix(c("albers_conical_equal_area", "aea", "azimuthal_equidistant", "aeqd", "lambert_cylindrical_equal_area", "cea", "lambert_azimuthal_equal_area", "laea", "lambert_conformal_conic", "lcc", "latitude_longitude", "longlat", "mercator", "merc", "orthographic", "ortho", "polar_stereographic", "stere", "stereographic", "stere", "transverse_mercator", "tmerc"), ncol=2, byrow=TRUE)

	m <- matrix(c("grid_mapping_name", "+proj", "false_easting", "+x_0","false_northing", "+y_0", "scale_factor_at_projection_origin", "+k_0", "scale_factor_at_central_meridian", "+k_0", "standard_parallel", "+lat_1", "standard_parallel1", "+lat_1", "standard_parallel2", "+lat_2", "longitude_of_central_meridian", "+lon_0", "longitude_of_projection_origin", "+lon_0", "latitude_of_projection_origin", "+lat_0", "straight_vertical_longitude_from_pole", "+lon_0",
	"longitude_of_prime_meridian", "+pm", "semi_major_axis", "+a", "semi_minor_axis", "+b", "inverse_flattening", "+rf",
	"earth_radius", "+a"), ncol=2, byrow=TRUE)

	# add logic that if prime merid is defined but not centr merid. centr merid is same as prime.

	i <- match(vars, m[,1])
	if (all(is.na(i))) {
		gg <- cbind(vars, vals)
		mtxt <- paste(apply(gg, 1, function(x) paste(x, collapse='=')), collapse='; ')
		warning("cannot process the crs\n", mtxt)
		return(NA)
	} else if (any(is.na(i))) {
		vr <- vars[is.na(i)]
		vl <- vals[is.na(i)]
		gg <- cbind(vr, vl)
		gg <- gg[!(gg[,1] %in% c("crs_wkt", "esri_pe_string")), ,drop=FALSE]
		if (NROW(gg) > 0) {
			mtxt <- paste(apply(gg, 1, function(x) paste(x, collapse='=')), collapse='\n')
			warning("cannot process these parts of the crs:\n", mtxt)
		}
		vars <- vars[!is.na(i)]
		vals <- vals[!is.na(i)]
		i <- stats::na.omit(i)
	}
	tab <- cbind(m[i,], vals)
	rr <- which(tab[,1] == "earth_radius")
	if (length(rr) > 0) {
		bb <- tab[rr,]
		bb[2] <- "+b"
		tab <- rbind(tab, bb)
	}
	p <- which(tab[,2] == '+proj')
	if (length(p) == 0) {
		warning("cannot create a valid crs\n", mtxt)
		return(NA)
	} else {
		tab <- rbind(tab[p, ], tab[-p, ])
	}
	j <- match(tab[1,3], prj[,1])
	tab[1,3] <- prj[j,2]
	cr <- paste(apply(tab[,2:3], 1, function(x) paste(x, collapse='=')), collapse=' ')
	crtst <- try(rast(crs=cr), silent=TRUE)
	if ( inherits(crtst, "try-error")) {
		mtxt <- paste(m, collapse='; ')
		warning("cannot create a valid crs\n", mtxt)
		return(NA)
	} else {
		return(cr)
	}
}


.ncdfTime <- function(nc, zvar, dim3, zval) {
	dodays <- TRUE
	dohours <- FALSE
	doseconds <- FALSE

	un <- nc$var[[zvar]]$dim[[dim3]]$units
	if (substr(un, 1, 10) == "days since") {
		startDate = as.Date(substr(un, 12, 22))
	} else if (substr(un, 1, 11) == "hours since") {
		dohours <- TRUE
		dodays <- FALSE
		startTime <- substr(un, 13, 30)
		mult <- 3600
	} else if (substr(un, 1, 13) == "seconds since") {
		doseconds <- TRUE
		dodays <- FALSE
		startTime = as.Date(substr(un, 15, 31))
		mult <- 1
	} else if (substr(un, 1, 12) == "seconds from") {
		doseconds <- TRUE
		dodays <- FALSE
		startTime = as.Date(substr(un, 14, 31))
		mult <- 1
	} else {
	  return(NULL)
	}
	if (!dodays) {
		start <- strptime(startTime, "%Y-%m-%d %H:%M:%OS", tz = "UTC")
		if (is.na(start)) start <- strptime(startTime, "%Y-%m-%d", tz = "UTC")
		if (is.na(start)) return(x)
		startTime <- start
		time <- startTime + as.numeric(zval) * mult
		time <- as.character(time)
		if (!is.na(time[1])) {
			return(time)
		}
	} else if (dodays) {
		# cal = nc$var[[zvar]]$dim[[dim3]]$calendar ?
		cal <- ncdf4::ncatt_get(nc, "time", "calendar")
		if (! cal$hasatt ) {
			greg <- TRUE
		} else {
			cal <- cal$value
			if (cal =='gregorian' | cal =='proleptic_gregorian' | cal=='standard') {
				greg <- TRUE
			} else if (cal == 'noleap' | cal == '365 day' | cal == '365_day') {
				greg <- FALSE
				nday <- 365
			} else if (cal == '360_day') {
				greg <- FALSE
				nday <- 360
			} else {
				greg <- TRUE
				warning('assuming a standard calender:', cal)
			}
		}
		if (greg) {
			time <- as.Date(time, origin=startDate)
		} else {
			startyear <-  as.numeric( format(startDate, "%Y") )
			startmonth <- as.numeric( format(startDate, "%m") )
			startday <- as.numeric( format(startDate, "%d") )
			year <- trunc( as.numeric(time)/nday )
			doy <- (time - (year * nday))
			origin <- paste(year+startyear, "-", startmonth, "-", startday, sep='')
			time <- as.Date(doy, origin=origin)
		}
		return(time)
	}
	return(NULL)
}



pointsCDF <- function(filename, varname, polygons=FALSE) {

	if (!("ncdf4" %in% rownames(utils::installed.packages()))) {
		warn("rast", "GDAL did not find an extent. installing the ncdf4 package may help")
		return(x)
	}

	zvar <- .varName(nc, varname, warn=TRUE)

	nc <- ncdf4::nc_open(filename, readunlim=FALSE, suppress_dimvals = TRUE)
	on.exit( ncdf4::nc_close(nc) )
	dims <- 1:3
	ncols <- nc$var[[zvar]]$dim[[dims[1]]]$len
	nrows <- nc$var[[zvar]]$dim[[dims[2]]]$len

	xx <- try(ncdf4::ncvar_get(nc, nc$var[[zvar]]$dim[[dims[1]]]$name), silent = TRUE)
	if (inherits(xx, "try-error")) {
		error("pointsCDF", "no x coordinates found")
	}
	yy <- try(ncdf4::ncvar_get(nc, nc$var[[zvar]]$dim[[dims[2]]]$name), silent = TRUE)
	if (inherits(yy, "try-error")) {
		error("pointsCDF", "no x coordinates found")
	}

	a <- ncdf4::ncatt_get(nc, zvar, "grid_mapping")
	prj <- NA
	if ( a$hasatt ) {
		try(atts <- ncdf4::ncatt_get(nc, a$value), silent=TRUE)
		try(prj <- .getCRSfromGridMap4(atts), silent=TRUE)
	}

	dim3 <- dims[3]
	dim3_vals <- try(ncdf4::ncvar_get(nc, nc$var[[zvar]]$dim[[dim3]]$name), silent = TRUE)
	if (inherits(dim3_vals, "try-error")) {
		dim3_vals <- seq_len(nc$var[[zvar]]$dim[[dim3]]$len)
	}
	nms <- NULL
	if ( nc$var[[zvar]]$dim[[dim3]]$name == "time" ) {
		try( nms <- .ncdfTime(nc, zvar, dim3, dim3_vals) )
	}

	d <- ncdf4::ncvar_get( nc, varid=zvar)
	nl <- dim(d)[3]
	v <- sapply(1:nl, function(i) d[,,i])

	natest1 <- ncdf4::ncatt_get(nc, zvar, "_FillValue")
	natest2 <- ncdf4::ncatt_get(nc, zvar, "missing_value")
	if (natest1$hasatt) {
		v[v==natest1$value] <- NA
	} else if (natest2$hasatt) {
		v[v==natest2$value] <- NA
	}
	if (!is.null(nms)) {
		colnames(v) <- nms
	}
	vect(cbind(rep(xx, length(yy)), rep(yy, each=length(xx))), atts=v, crs=prj)
}

