
.ncdf_extent <- function(x) {

	stopifnot(requireNamespace("ncdf4"))
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

