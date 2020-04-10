

gdal_info <- function(filename, print=TRUE, opt1="", opt2="") {
	x <- .gdalinfo(f, opt1, opt2)
	if (print) {
		cat(x)
		invisible(x)
	} else {
		return(x)
	}
}

